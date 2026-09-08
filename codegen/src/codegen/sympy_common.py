from collections import namedtuple
from collections.abc import Callable, Sequence

import numpy as np

import sympy as sym
from sympy import Basic, Symbol, Matrix, ImmutableMatrix, MatrixSymbol
from sympy.utilities.iterables import numbered_symbols
from sympy.codegen.ast import Assignment
from sympy.printing.cxx import CXX17CodePrinter

NamedExpr = namedtuple("NamedExpr", ["name", "expr"])


def make_vector(name: str, dim: int, **kwargs) -> Matrix:
    return Matrix([[Symbol(f"{name}[{i}]", **kwargs)] for i in range(dim)])


def make_matrix(name: str, rows: int, cols: int, **kwargs) -> Matrix:
    return Matrix(
        [
            [Symbol(f"{name}[{i},{j}]", **kwargs) for j in range(cols)]
            for i in range(rows)
        ]
    )


def name_expr(name: str, expr: Basic) -> NamedExpr:
    if hasattr(expr, "shape"):
        s = sym.MatrixSymbol(name, *expr.shape)
    else:
        s = Symbol(name)
    return NamedExpr(s, expr)


def find_by_name(name_exprs: Sequence[NamedExpr], name: str) -> NamedExpr | None:
    return next(
        (name_expr for name_expr in name_exprs if str(name_expr[0]) == name), None
    )


def explicit(x: Basic) -> Basic:
    """Explicit matrix view of a MatrixSymbol, passing anything else through."""
    return x.as_explicit() if hasattr(x, "as_explicit") else x


def _collapse(x: Basic) -> Basic:
    """Force a matrix expression down to its entries.

    Substituting a `MatrixSymbol` by an explicit matrix leaves an unevaluated
    `MatAdd`/`MatMul` behind, which the next substitution would have to walk
    again.
    """
    if hasattr(x, "doit"):
        x = x.doit()
    if hasattr(x, "as_explicit"):
        x = x.as_explicit()
    return ImmutableMatrix(x) if hasattr(x, "shape") else x


class Derivation:
    """A straight-line sequence of named intermediate expressions.

    Each is in terms of the ones before it, which is what the printer turns
    into C++ locals. Keeping them here rather than in a bare list buys
    ``resolve`` and ``check_same``.
    """

    def __init__(self) -> None:
        self.name_exprs: list[NamedExpr] = []
        self._by_name: dict[Basic, Basic] = {}
        self._at_cache: tuple[int | None, dict[Basic, Basic]] = (None, {})

    def add(self, name: str, expr: Basic) -> NamedExpr:
        """Name an expression and append it to the derivation.

        @param name is the C++ identifier the expression will be emitted as
        @return the NamedExpr, whose ``.name`` refers to it downstream
        """
        if isinstance(expr, Matrix):
            expr = ImmutableMatrix(expr)
        ne = name_expr(name, expr)
        self.name_exprs.append(ne)
        self._by_name[ne.name] = ne.expr
        return ne

    def record(self, named_expr: NamedExpr) -> NamedExpr:
        """Adopt an already-built NamedExpr, so `resolve` can see through it.

        @return the same pair, for chaining
        """
        self.name_exprs.append(named_expr)
        self._by_name[named_expr.name] = named_expr.expr
        return named_expr

    def resolve(self, expr: Basic, at: dict | None = None) -> Basic:
        """Substitute every named intermediate away, down to the inputs.

        Backwards without @p at, which the dependency order makes a single
        pass; forwards with it, so each definition collapses to a number as it
        is reached and no expression grows past one definition.

        @param at optionally evaluates every definition at that point
        @return the closed form of @p expr, or its value at @p at
        """
        if at is not None:
            # the same point is used for every resolve of a derivation, so
            # evaluate each definition once rather than per call
            key = id(at)
            if self._at_cache[0] != key:
                values: dict[Basic, Basic] = {}
                for ne in self.name_exprs:
                    values[ne.name] = _collapse(ne.expr.subs(at).subs(values))
                self._at_cache = (key, values)
            expr = _collapse(expr.subs(at).subs(self._at_cache[1]))
        else:
            for i in range(len(self.name_exprs) - 1, -1, -1):
                ne = self.name_exprs[i]
                expr = expr.subs(ne.name, ne.expr)
        if any(expr.has(ne.name) for ne in self.name_exprs):
            raise RuntimeError("named expressions did not resolve")
        return expr

    def check_same(self, what: str, expr_a: Basic, expr_b: Basic) -> None:
        """Fail unless two expressions agree once fully resolved.

        @param what names the shortcut, for the error message
        """
        diff = sym.simplify(sym.expand(self.resolve(expr_a) - self.resolve(expr_b)))
        if any(e != 0 for e in diff):
            raise AssertionError(f"{what}: shortcut does not match chain rule\n{diff}")


class StructuredMatrix:
    """A symbolic matrix declared as a grid of entry classes.

    A number is a structural constant, a string names the entry's class and
    becomes a symbol.  The index sets the classes define are read back off the
    grid rather than restated.
    """

    def __init__(self, name, spec):
        cols = {len(row) for row in spec}
        if len(cols) != 1:
            raise ValueError(f"{name}: rows of unequal length")

        self.rows = len(spec)
        self.cols = cols.pop()
        self.size = self.rows * self.cols
        self._spec = [list(row) for row in spec]

        symbols = MatrixSymbol(name, self.rows, self.cols)
        self.matrix = Matrix(
            [
                [
                    symbols[i, j] if isinstance(entry, str) else entry
                    for j, entry in enumerate(row)
                ]
                for i, row in enumerate(self._spec)
            ]
        )

    def entries(self, *classes):
        """(row, column) of every entry in `classes`, in storage order."""
        return [
            (i, j)
            for j in range(self.cols)
            for i in range(self.rows)
            if self._spec[i][j] in classes
        ]

    def flat_index(self, row, col):
        """Index of an entry in the column major storage."""
        return row + self.rows * col


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

    def _as_ordered_terms(self, expr, order=None):
        # A compiler may only contract a multiply into a neighbouring add
        # within one expression, and it associates left to right, so a sum
        # folds into fused multiply-adds only if it does not *start* with a
        # product.  sympy's canonical order routinely puts one first:
        #
        #   -H1[1]*Tp2[2] + H1[2]*Tp2[1] + M[20]   ->  fmul, fmsub, fadd
        #
        # while the same sum led by its plain term costs one operation less:
        #
        #   M[20] - H1[1]*Tp2[2] + H1[2]*Tp2[1]   ->  fmsub, fmadd
        #
        # Hoisting the non-product terms to the front is enough to get the
        # second form.  It changes the order of floating point additions, so
        # results move in the last bits -- which is why it is done here, in
        # the printer, and not by rewriting the expressions themselves.
        terms = super()._as_ordered_terms(expr, order=order)
        plain = [t for t in terms if not t.is_Mul]
        if not plain:
            return terms
        return plain + [t for t in terms if t.is_Mul]

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
    printer: CXX17CodePrinter,
    name_exprs: Sequence[NamedExpr],
    outputs: Sequence[Basic],
    run_cse: bool = True,
    pre_expr_hook: Callable[[Basic], str | None] | None = None,
    post_expr_hook: Callable[[Basic], str | None] | None = None,
) -> str:
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
                # std::array: same layout and the same `x[i]` access, but it
                # carries its own size and cannot decay to a pointer
                lines.append(f"std::array<T, {np.prod(expr.shape)}> {var}{{}};")
                lines.extend(code.split("\n"))
            else:
                lines.append("const auto " + code)
        else:
            # A scalar output is a reference parameter, so it is assigned
            # directly; only a matrix output is written through its span.
            lines.extend(code.split("\n"))

        if post_expr_hook is not None:
            code = post_expr_hook(var)
            if code is not None:
                lines.extend(code.split("\n"))

    return "\n".join(lines)


def my_function_print(
    printer: CXX17CodePrinter,
    name: str,
    inputs: Sequence[Basic],
    name_exprs: Sequence[NamedExpr],
    outputs: Sequence[Basic],
    run_cse: bool = True,
) -> str:
    """Emit a whole function, signature included, in the ACTS dialect.

    Matrices are passed as sized spans and scalars by reference, so a caller
    cannot hand over the wrong extent. The detray backend spells its own
    signatures in ``gen_cxx_code``.
    """

    def extent(sym):
        return int(np.prod(sym.shape))

    def input_param(input):
        if isinstance(input, MatrixSymbol):
            return f"std::span<const T, {extent(input)}> {input.name}"
        return f"const T {input.name}"

    def output_param(output):
        if isinstance(output, MatrixSymbol):
            return f"std::span<T, {extent(output)}> {output.name}"
        return f"T& {output.name}"

    lines = []

    params = [input_param(input) for input in inputs] + [
        output_param(output) for output in outputs
    ]
    head = f"template<typename T> void {name}({', '.join(params)}) {{"

    lines.append(head)

    code = my_expression_print(printer, name_exprs, outputs, run_cse=run_cse)
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)
