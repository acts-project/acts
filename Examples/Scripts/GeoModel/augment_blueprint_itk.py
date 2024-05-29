import sqlite3
import os
import argparse


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-i",
        "--input",
        type=str,
        default="",
        help="Input original SQLite database",
    )

    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Output SQLite database",
    )

    args = p.parse_args()

    # Copy input to output
    cmd = f"cp {args.input} {args.output}"
    print(">> Performing", cmd)
    os.system(cmd)

    # Creating the database connection
    connection = sqlite3.connect(args.output)
    cursor = connection.cursor()

    print(">> Creating new table called Blueprint")

    cursor.execute(
        "CREATE TABLE Blueprint(id INT, type TEXT, name TEXT, bounds TEXT, internals TEXT, binnings TEXT, materials TEXT)"
    )

    # Bounds nomenclature is:
    # (type;[values])
    # type: cyl, box
    # values:
    #  - cyl: [rmin, rmax, zmin, zmax]
    #  - box: [xmin, xmax, ymin, ymax, zmin, zmax]
    # special characters: i ... inherited (from mother),
    #                     c ... calculated (from container siblings),
    #                     e ... envelope (from internals) with default value
    #                     eX ... envelope (from internals) with value X
    #                     * ... any value
    #
    # Binning nomenclature is:

    print(">> Filling the entries ...")

    # Augmenting the GeoModel sqlite
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (0, 'root', 'ITk', 'cyl;0.,1180,-3500.,3500.','children;ITk/NegSector,ITk/Central,ITk/PosSector','z', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (100, 'leaf', 'ITk/NegSector', 'cyl;i,i,-3500.,-3050.','','', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (101, 'container', 'ITk/Central', 'cyl;i,i,-3050.,3050.',
            'children;ITk/Central/Detectors,ITk/Central/Outer','r', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (102, 'leaf', 'ITk/PosSector', 'cyl;,i,i,3050.,3500','','', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (91, 'container', 'ITk/Central/Detectors', 'cyl;i,1070,i,i',
            'children;ITk/Central/Detectors/Pixels,ITk/Central/Detectors/Strips','r', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (92, 'leaf', 'ITk/Central/Outer', 'cyl;1070,i,i,i','','', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (81, 'container', 'ITk/Central/Detectors/Pixels', 'cyl;i,345,i,i',
            'children;ITk/Central/Detectors/Pixels/InnerPixels,ITk/Central/Detectors/Pixels/OuterPixels','r', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (82, 'container', 'ITk/Central/Detectors/Strips', 'cyl;345.,i,i,i',
            'children;ITk/Central/Detectors/Strips/NegSector,ITk/Central/Detectors/Strips/Barrel,ITk/Central/Detectors/Strips/PosSector',
            'z', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (61, 'leaf', 'ITk/Central/Detectors/Pixels/InnerPixels', 'cyl;i,130,i,i','','', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (62, 'leaf', 'ITk/Central/Detectors/Pixels/OuterPixels', 'cyl;130,i,i,i','','', '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (71, 'container', 'ITk/Central/Detectors/Strips/NegSector', 'cyl;i,i,i,-1400.',
            'children;ITk/Central/Detectors/Strips/NegSector/NegEndcap,ITk/Central/Detectors/Strips/NegSector/OuterGap',
            'r', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (62, 'leaf', 'ITk/Central/Detectors/Strips/Barrel', 'cyl;i,i,-1400.,1400.','','', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (73, 'container', 'ITk/Central/Detectors/Strips/PosSector', 'cyl;i,i,1400.,i',
            'children;ITk/Central/Detectors/Strips/PosSector/PosEndcap,ITk/Central/Detectors/Strips/PosSector/OuterGap',
            'r', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (61, 'leaf', 'ITk/Central/Detectors/Strips/NegSector/NegEndcap', 'cyl;i,990.,i,i',
            '',
            '', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (72, 'leaf', 'ITk/Central/Detectors/Strips/NegSector/OuterGap', 'cyl;990.,i,i,i',
            '',
            '', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (63, 'leaf', 'ITk/Central/Detectors/Strips/PosSector/PosEndcap', 'cyl;i,990,i,i',
            '',
            '', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (74, 'leaf', 'ITk/Central/Detectors/Strips/PosSector/OuterGap', 'cyl;990.,i,i,i',
            '',
            '', '')
    """
    )

    connection.commit()

    print(">> Commited and done.")
