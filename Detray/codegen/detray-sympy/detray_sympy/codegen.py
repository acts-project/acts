from detray_sympy.common import (
    cxx_printer,
    my_expression_print,
)


def render_dim_requirement(name, j):
    if not hasattr(j, "shape"):
        return "(detray::concepts::scalar<%s_t>)" % (name)
    # HACK: This is for the 8x1 path-to-free matrix.
    elif len(j.shape) == 1 or (
        len(j.shape) == 2 and j.shape[1] == 1 and j.shape[0] <= 3
    ):
        if j.shape[0] == 3:
            return "(detray::concepts::vector3D<%s_t>)" % (name)
        else:
            return (
                "(detray::concepts::vector<%s_t> && detray::traits::rows<%s_t> == %d)"
                % (name, name, j.shape[0])
            )
    elif len(j.shape) == 2:
        if j.shape[0] == j.shape[1]:
            return (
                "(detray::concepts::square_matrix<%s_t> && detray::traits::max_rank<%s_t> == %d)"
                % (name, name, j.shape[0])
            )
        elif j.shape[0] == 1:
            return (
                "(detray::concepts::row_matrix<%s_t> && detray::traits::columns<%s_t> == %d)"
                % (name, name, j.shape[1])
            )
        else:
            return (
                "(detray::concepts::matrix<%s_t> && detray::traits::rows<%s_t> == %d && detray::traits::columns<%s_t> == %d)"
                % (name, name, j.shape[0], name, j.shape[1])
            )


def gen_cxx_code(function_name, inputs, outputs, run_cse=True, printer=None):
    if printer is None:
        printer = cxx_printer

    template_types = []

    for i, j in inputs:
        template_types.append("%s_t" % i)
    for i, j in outputs:
        template_types.append("%s_t" % i)

    lines = []

    lines.append(
        "template <%s>" % (", ".join("typename %s" % s for s in template_types))
    )
    lines.append("DETRAY_HOST_DEVICE void inline %s (" % function_name)
    lines.append(", ".join("const %s_t & %s" % (i, i) for i, _ in inputs) + ",")
    lines.append(", ".join("%s_t & %s" % (i, i) for i, _ in outputs))
    lines.append(")")
    if len(inputs) > 0 or len(outputs) > 0:
        lines.append("requires(")
        lines.append(
            " && ".join(render_dim_requirement(i, j) for i, j in inputs + outputs)
        )
        lines.append(")")
    lines.append("{")

    for i, j in inputs:
        if not hasattr(j, "shape"):
            continue
        # HACK: This is for the 8x1 path-to-free matrix.
        if len(j.shape) == 1 or (
            len(j.shape) == 2 and j.shape[1] == 1 and j.shape[0] <= 3
        ):
            for k in range(j.shape[0]):
                if j[k] == 0:
                    lines.append(
                        "assert((getter::element<{idx}>({var}) == 0.f));".format(
                            var=i, idx=k
                        )
                    )
                elif j[k] == 1:
                    lines.append(
                        "assert((getter::element<{idx}>({var}) == 1.f));".format(
                            var=i, idx=k
                        )
                    )
        elif len(j.shape) == 2:
            for k in range(j.shape[0]):
                for l in range(j.shape[1]):
                    if j[k, l] == 0:
                        lines.append(
                            "assert((getter::element<{idx1}, {idx2}>({var}) == 0.f));".format(
                                var=i, idx1=k, idx2=l
                            )
                        )
                    elif j[k, l] == 1:
                        lines.append(
                            "assert((getter::element<{idx1}, {idx2}>({var}) == 1.f));".format(
                                var=i, idx1=k, idx2=l
                            )
                        )

    code = my_expression_print(
        printer,
        outputs,
        [x[0] for x in outputs],
        run_cse=run_cse,
    )
    lines.extend([f"  {l}" for l in code.split("\n")])

    lines.append("}")

    return "\n".join(lines)
