from collections import namedtuple

import numpy as np

import sympy as sym
from sympy import Symbol, Matrix, ImmutableMatrix, MatrixSymbol
from sympy.utilities.iterables import numbered_symbols
from sympy.codegen.ast import Assignment
from sympy.printing.cxx import CXX17CodePrinter

NamedExpr = namedtuple("NamedExpr", ["name", "expr"])


def make_vector(name, dim, **kwargs):
    return Matrix([[Symbol(f"{name}[{i}]", **kwargs)] for i in range(dim)])


def make_matrix(name, rows, cols, **kwargs):
    return Matrix(
        [
            [Symbol(f"{name}[{i},{j}]", **kwargs) for j in range(cols)]
            for i in range(rows)
        ]
    )


def name_expr(name, expr):
    if hasattr(expr, "shape"):
        s = sym.MatrixSymbol(name, *expr.shape)
    else:
        s = Symbol(name)
    return NamedExpr(s, expr)


def find_by_name(name_exprs, name):
    return next(
        (name_expr for name_expr in name_exprs if str(name_expr[0]) == name), None
    )


# Matrix-element access strategies. A strategy is a callable
# ``(printer, MatrixElement) -> str`` that renders how a single matrix entry is
# accessed in the generated C++. This is the one genuine point of divergence
# between the ACTS and detray backends, so it is made pluggable here while the
# rest of the printer (and all of the symbolic machinery below) stays shared.


def flat_index_accessor(printer, expr):
    """ACTS dialect: flat column-major buffer indexing ``var[i + j*rows]``."""
    from sympy.printing.precedence import PRECEDENCE

    return "{}[{}]".format(
        printer.parenthesize(expr.parent, PRECEDENCE["Atom"], strict=True),
        expr.i + expr.j * expr.parent.shape[0],
    )


def getter_element_accessor(printer, expr):
    """detray dialect: algebra-abstracted ``getter::element<i, j>(var)``."""
    from sympy.printing.precedence import PRECEDENCE

    var = printer.parenthesize(expr.parent, PRECEDENCE["Atom"], strict=True)
    # HACK: This is for the 8x1 path-to-free matrix.
    if expr.parent.shape[1] == 1 and expr.parent.shape[0] <= 3:
        return "getter::element<{idx}>({var})".format(var=var, idx=expr.i)
    return "getter::element<{idx1}, {idx2}>({var})".format(
        var=var, idx1=expr.i, idx2=expr.j
    )


class MyCXXCodePrinter(CXX17CodePrinter):
    def __init__(self, *args, element_accessor=None, **kwargs):
        super().__init__(*args, **kwargs)
        # Default to the ACTS flat-index dialect so existing ACTS generators
        # are byte-for-byte unaffected.
        self._element_accessor = element_accessor or flat_index_accessor

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for j in range(cols) for i in range(rows))

    def _print_MatrixElement(self, expr):
        return self._element_accessor(self, expr)

    def _print_Pow(self, expr):
        from sympy.core.numbers import equal_valued, Float
        from sympy.codegen.ast import real

        suffix = self._get_func_suffix(real)
        if equal_valued(expr.exp, -0.5):
            return "%s/%ssqrt%s(%s)" % (
                self._print_Float(Float(1.0)),
                self._ns,
                suffix,
                self._print(expr.base),
            )
        return super()._print_Pow(expr)


class MyCXXCodePrinterWithoutKnownAssignment(MyCXXCodePrinter):
    def _print_Assignment(self, expr):
        if expr.rhs == 0 or expr.rhs == 1:
            return ""
        return super()._print_Assignment(expr)


cxx_printer = MyCXXCodePrinter()
cxx_printer_wo_known = MyCXXCodePrinterWithoutKnownAssignment()


def inflate_expr(name_expr):
    name, expr = name_expr

    result = []
    references = []

    if hasattr(expr, "shape"):
        for indices in np.ndindex(expr.shape):
            result.append((name[indices], expr[indices]))
            references.append((name, expr.shape, indices))
    else:
        result.append((name, expr))
        references.append(None)

    return result, references


def inflate_exprs(name_exprs):
    result = []
    references = []
    for name_expr in name_exprs:
        res, refs = inflate_expr(name_expr)
        result.extend(res)
        references.extend(refs)
    return result, references


def deflate_exprs(name_exprs, references):
    result = []
    deflated = {}

    for name_expr, reference in zip(name_exprs, references):
        if reference is None:
            result.append(name_expr)
        else:
            _, expr = name_expr
            name, shape, indices = reference
            if name not in deflated:
                e = Matrix(np.zeros(shape))
                result.append(NamedExpr(name, e))
                deflated[name] = e
            deflated[name][*indices] = expr

    another_result = []
    for name_expr in result:
        name, expr = name_expr
        if isinstance(expr, Matrix):
            another_result.append(NamedExpr(name, ImmutableMatrix(expr)))
        else:
            another_result.append(name_expr)

    return another_result


def my_subs(expr, sub_name_exprs):
    sub_name_exprs, _ = inflate_exprs(sub_name_exprs)

    result = expr.expand()
    result = result.subs([(e, n) for n, e in sub_name_exprs])
    result = sym.simplify(result)
    return result


def build_dependency_graph(name_exprs):
    graph = {}
    for name, expr in name_exprs:
        graph[name] = expr.free_symbols
    return graph


def build_influence_graph(name_exprs):
    graph = {}
    for name, expr in name_exprs:
        for s in expr.free_symbols:
            graph.setdefault(s, set()).add(name)
    return graph


def order_exprs_by_input(name_exprs):
    all_expr_names = set().union(name for name, _ in name_exprs)
    all_expr_symbols = set().union(*[expr.free_symbols for _, expr in name_exprs])
    inputs = all_expr_symbols - all_expr_names

    order = {}

    order.update({i: 0 for i in inputs})

    while len(order) < len(inputs) + len(name_exprs):
        for name, expr in name_exprs:
            symbols_order = [order.get(s, None) for s in expr.free_symbols]
            if None in symbols_order:
                continue
            if len(symbols_order) == 0:
                order[name] = 0
            else:
                order[name] = max(symbols_order) + 1

    result = name_exprs
    result = sorted(result, key=lambda n_e: len(n_e[1].args))
    result = sorted(result, key=lambda n_e: len(n_e[1].free_symbols))
    result = sorted(result, key=lambda n_e: order[n_e[0]])
    return result


def order_exprs_by_output(name_exprs, outputs):
    name_expr_by_name = {name_expr[0]: name_expr for name_expr in name_exprs}

    def get_inputs(output):
        name_expr = name_expr_by_name.get(output, None)
        if name_expr is None:
            return set()
        inputs = set(name_expr[1].free_symbols)
        inputs.update(*[get_inputs(name) for name in inputs])
        return inputs

    result = []
    done = set()

    for output in outputs:
        inputs = get_inputs(output) - done
        result.extend(
            order_exprs_by_input(
                [name_exprs for name_exprs in name_exprs if name_exprs[0] in inputs]
            )
        )
        result.append(name_expr_by_name[output])
        done.update(inputs)
        done.add(output)

    return result


def my_cse(name_exprs, inflate_deflate=True, simplify=True):
    sub_symbols = numbered_symbols()

    if inflate_deflate:
        name_exprs, references = inflate_exprs(name_exprs)

    names = [x[0] for x in name_exprs]
    exprs = [x[1] for x in name_exprs]

    sub_exprs, simp_exprs = sym.cse(exprs, symbols=sub_symbols)

    if simplify:
        sub_exprs = [(n, sym.simplify(e)) for n, e in sub_exprs]
        simp_exprs = [sym.simplify(e) for e in simp_exprs]

    simp_name_exprs = list(zip(names, simp_exprs))
    if inflate_deflate:
        simp_name_exprs = deflate_exprs(simp_name_exprs, references)

    name_exprs = []
    name_exprs.extend(sub_exprs)
    name_exprs.extend(simp_name_exprs)

    return name_exprs


def my_expression_print(
    printer, name_exprs, outputs, run_cse=True, pre_expr_hook=None, post_expr_hook=None
):
    if run_cse:
        name_exprs = my_cse(name_exprs, inflate_deflate=True)
    name_exprs = order_exprs_by_output(name_exprs, outputs)

    lines = []

    for var, expr in name_exprs:
        if pre_expr_hook is not None:
            code = pre_expr_hook(var)
            if code is not None:
                lines.extend(code.split("\n"))

        code = printer.doprint(Assignment(var, expr))
        if var not in outputs:
            if hasattr(expr, "shape"):
                lines.append(f"T {var}[{np.prod(expr.shape)}];")
                lines.extend(code.split("\n"))
            else:
                lines.append("const auto " + code)
        else:
            if hasattr(expr, "shape"):
                lines.extend(code.split("\n"))
            else:
                lines.append("*" + code)

        if post_expr_hook is not None:
            code = post_expr_hook(var)
            if code is not None:
                lines.extend(code.split("\n"))

    return "\n".join(lines)


def my_function_print(printer, name, inputs, name_exprs, outputs, run_cse=True):
    def input_param(input):
        if isinstance(input, MatrixSymbol):
            return f"const T* {input.name}"
        return f"const T {input.name}"

    def output_param(name):
        return f"T* {name}"

    lines = []

    params = [input_param(input) for input in inputs] + [
        output_param(output) for output in outputs
    ]
    head = f"template<typename T> void {name}({", ".join(params)}) {{"

    lines.append(head)

    code = my_expression_print(printer, name_exprs, outputs, run_cse=run_cse)
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)
