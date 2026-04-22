import sys
import sympy
from sympy.tensor.array.expressions import ArraySymbol

from detray_sympy.common import name_expr
from detray_sympy.output import write_out_file
from detray_sympy.codegen import gen_cxx_code
import detray_sympy.matrices


def get_linear_transport_jacobian_matrix(gradient):
    matrix = detray_sympy.matrices.get_generic_matrix_D(gradient=gradient)

    num_unknown = 0

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i, j] != 0 and matrix[i, j] != 1:
                num_unknown += 1

    array = ArraySymbol("A", (num_unknown,))

    i_unknown = 0

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i, j] != 0 and matrix[i, j] != 1:
                matrix[i, j] = array[i_unknown]
                i_unknown += 1

    return array, matrix


def gen_code(gradient):
    array, matrix = get_linear_transport_jacobian_matrix(gradient)

    lines = []

    struct_name = (
        "transport_jacobian_matrix_with_gradient"
        if gradient
        else "transport_jacobian_matrix_without_gradient"
    )

    lines.append("template <concepts::algebra algebra_t>")

    lines.append("struct %s {" % struct_name)

    lines.append("using algebra_type = algebra_t;")
    lines.append("using scalar_type = dscalar<algebra_type>;")
    lines.append("using value_type = scalar_type;")

    lines.append("template<std::size_t I, std::size_t J>")
    lines.append(
        "DETRAY_HOST_DEVICE scalar_type element() const requires (I < 8 && J < 8) {"
    )

    first = True
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            stmt = "if" if first else "else if"
            first = False

            lines.append("%s constexpr (I == %d && J == %d) {" % (stmt, i, j))

            if matrix[i, j] == 0:
                lines.append("return static_cast<scalar_type>(0.f);")
            elif matrix[i, j] == 1:
                lines.append("return static_cast<scalar_type>(1.f);")
            else:
                lines.append("return m_contents[%d];" % (matrix[i, j].indices[0]))

            lines.append("}")

    lines.append("}")

    lines.append("template<std::size_t I, std::size_t J>")
    lines.append("DETRAY_HOST_DEVICE scalar_type& element() requires(")

    lines.append(
        " || ".join(
            "(I == %d && J == %d)" % (i, j)
            for i in range(matrix.shape[0])
            for j in range(matrix.shape[1])
            if (matrix[i, j] != 0 and matrix[i, j] != 1)
        )
    )

    lines.append(") {")

    first = True
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i, j] == 0 or matrix[i, j] == 1:
                continue
            stmt = "if" if first else "else if"
            first = False

            lines.append(
                "%s constexpr (I == %d && J == %d) { return m_contents[%d]; }"
                % (stmt, i, j, matrix[i, j].indices[0])
            )

    lines.append("}")

    lines.append(
        "DETRAY_HOST_DEVICE explicit operator dmatrix<algebra_t, 8, 8>() const {"
    )
    lines.append("dmatrix<algebra_t, 8, 8> rv;")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            lines.append(
                "getter::element<%d, %d>(rv) = element<%d, %d>();" % (i, j, i, j)
            )

    lines.append("return rv;")
    lines.append("}")

    lines.append("DETRAY_HOST_DEVICE static constexpr %s identity() {" % struct_name)
    lines.append("%s rv;" % struct_name)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if i == j:
                assert matrix[i, j] != 0
                if matrix[i, j] != 1:
                    lines.append(
                        "rv.element<%d, %d>() = static_cast<scalar_type>(1.f);" % (i, j)
                    )
            else:
                assert matrix[i, j] != 1
                if matrix[i, j] != 0:
                    lines.append(
                        "rv.element<%d, %d>() = static_cast<scalar_type>(0.f);" % (i, j)
                    )

    lines.append("return rv;")
    lines.append("}")

    lines.append("private:")
    lines.append("std::array<scalar_type, %d> m_contents;" % array.shape[0])

    lines.append("};\n")

    return "\n".join(lines)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        output = sys.argv[1]
    else:
        output = None

    c1 = gen_code(gradient=False)
    c2 = gen_code(gradient=True)

    write_out_file(c1 + c2, output)
