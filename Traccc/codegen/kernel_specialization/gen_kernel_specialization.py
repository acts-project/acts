import string
import argparse
import pathlib

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "template",
        type=pathlib.Path,
        help="the template path",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        help="the output path to write to",
        required=True,
    )

    parser.add_argument(
        "--detector",
        type=str,
        help="the name of the detector",
    )

    parser.add_argument(
        "--bfield",
        type=str,
        help="the name of the bfield",
    )

    parser.add_argument(
        "--model",
        type=str,
        help="the name of the programming mode",
        default="cpu",
    )

    args = parser.parse_args()

    with open(args.template, "r") as f:
        src = string.Template(f.read())

    subs = {}

    subs["SOURCE_DIR"] = args.template.parent

    det = getattr(args, "detector", None)
    if det is not None:
        subs["DETECTOR_NAME"] = det

    bfield = getattr(args, "bfield", None)
    if bfield is not None:
        # HACK: Perform some transformations to make the types line up...
        # This could be resolved in the C++ file itself in the future.
        bfield_name = bfield

        if args.model != "cpu" and bfield != "const":
            bfield_name = args.model + "::" + bfield_name

        bfield_name += "_bfield_backend_t"

        if bfield != "inhom_texture":
            bfield_name += "<scalar>"

        subs["BFIELD_NAME"] = bfield_name

    result = src.substitute(subs)

    with open(args.output, "w") as f:
        f.write(result)
