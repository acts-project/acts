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
            (37, 
            'root', 
            'ITk', 
            'cyl;0.,1180,-3500.,3500.',
            'children;NegSector,Central,PosSector',
            'z', 
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (36, 
            'leaf', 
            'ITk/NegSector', 
            'cyl;i,i,-3500.,-3050.',
            '',
            '',
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (35, 
            'container', 
            'ITk/Central', 
            'cyl;i,i,-3050.,3050.',
            'children;BeamPipe,Detectors,Outer',
            'r', 
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (33, 
            'leaf', 
            'ITk/PosSector',
            'cyl;,i,i,3050.,3500',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (0, 
            'leaf',
            'ITk/Central/BeamPipe',
            'cyl;i,25,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (32, 
            'container',
            'ITk/Central/Detectors',
            'cyl;25,1070,i,i',
            'children;Pixels,Strips',
            'r', 
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (34, 
            'leaf', 
            'ITk/Central/Outer', 
            'cyl;1070,i,i,i',
            '',
            '',
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (30, 
            'container',
            'ITk/Central/Detectors/Pixels',
            'cyl;i,345,i,i',
            'children;InnerPixels,OuterPixels',
            'r',
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (31, 
            'container',
            'ITk/Central/Detectors/Strips',
            'cyl;345.,i,i,i',
            'children;NegSector,Barrel,PosSector',
            'z', 
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (20, 
            'container',
            'ITk/Central/Detectors/Pixels/InnerPixels',
            'cyl;i,130,i,i',
            'children;NegEndcap,Barrel,PosEndcap',
            'z',
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (25, 
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels',
            'cyl;130,i,i,i',
            'children;NegEndcap,NegInclined,Barrel,PosInclined,PosEndcap',
            'z',
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (27,
            'container',
            'ITk/Central/Detectors/Strips/NegSector',
            'cyl;i,i,i,-1400.',
            'children;NegEndcap,OuterGap',
            'r', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (18,
            'leaf',
            'ITk/Central/Detectors/Strips/Barrel',
            'cyl;i,i,-1400.,1400.',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (29,
            'container',
            'ITk/Central/Detectors/Strips/PosSector',
            'cyl;i,i,1400.,i',
            'children;PosEndcap,OuterGap',
            'r', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (17,
            'leaf',
            'ITk/Central/Detectors/Strips/NegSector/NegEndcap',
            'cyl;i,990.,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (26,
            'leaf',
            'ITk/Central/Detectors/Strips/NegSector/OuterGap',
            'cyl;990.,i,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (19,
            'leaf',
            'ITk/Central/Detectors/Strips/PosSector/PosEndcap', 
            'cyl;i,990,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (28,
            'leaf',
            'ITk/Central/Detectors/Strips/PosSector/OuterGap',
            'cyl;990.,i,i,i',
            '',
            '', '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (21, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap', 
            'cyl;i,i,i,-1100.',
            'children;Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (22, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined',
            'cyl;i,i,-1100.,-380.',
            'children;Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (10, 
            'leaf', 
            'ITk/Central/Detectors/Pixels/OuterPixels/Barrel',
            'cyl;i,i,-380.,380.',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (23, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined',
            'cyl;i,i,380.,1100',
            'children;Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (24, 
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap',
            'cyl;i,i,1100.,i',
            'children;Ring0,Ring1,Ring2',
            'r',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (4,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring0', 'cyl;i,205.,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (5,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring1',
            'cyl;205.,255.,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (6, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring2',
            'cyl;255.,i,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (7, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring0',
            'cyl;i,205.,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (8,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring1',
            'cyl;205.,255.,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (9,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring2',
            'cyl;255.,i,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (11,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring0',
            'cyl;i,205.,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (12,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring1',
            'cyl;205.,255.,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (13, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring2',
            'cyl;255.,i,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (14, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring0',
            'cyl;i,205.,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (15,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring1',
            'cyl;205.,255.,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (16,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring2',
            'cyl;255.,i,i,i',
            '',
            '',
            '')
    """
    )
    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (1, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/NegEndcap',
            'cyl;i,i,i,-255.',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (2,
            'container',
            'ITk/Central/Detectors/Pixels/InnerPixels/Barrel',
            'cyl;i,i,-255.,255.',
            'children;Layer0,Layer1,Layer2,Layer3',
            'r', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (200,
            'container',
            'ITk/Central/Detectors/Pixels/InnerPixels/Barrel',
            'cyl;i,i,-255.,255.',
            'children;Layer0,*,Layer1,*',
            'r', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (201,
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/Barrel/Layer0',
            'cyl;i,30,i,i',
            'sensitives;kdt;i,30,i,i',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (202,
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/Barrel/Layer1',
            'cyl;80,120,i,i',
            'sensitives;kdt;80,120,i,i',
            '', 
            '')
    """
    )

    cursor.execute(
        """
    INSERT INTO Blueprint VALUES
            (3, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/PosEndcap',
            'cyl;i,i,255.,i',
            '',
            '',
            '')
    """
    )

    connection.commit()

    print(">> Commited and done.")
