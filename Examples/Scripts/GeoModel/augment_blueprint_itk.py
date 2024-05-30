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

    # --------------------------------------------------------------------------------------
    # ITk overall
    ITk_r_max = 1180.
    ITk_z_max = 3500.

    # ITK central
    ITK_central_z_max = 3050.
    ITk_central_r_max = 1070.

    # --------------------------------------------------------------------------------------
    # Beam Pipe section - BP
    BP_r_max = 27.

    # --------------------------------------------------------------------------------------
    # Inner Pixels sections - IP
    IP_r_min = BP_r_max
    IP_r_max = 124.
    IP_z_mid = 245.

    IP_b_z_max = 255.
   
    IP_coupled_rings_bins_r = 2 
    IP_coupled_rings_bins_phi = 20 

    IP_inner_rings_bins_r = 1 
    IP_inner_rings_bins_phi = 30 

    IP_outer_rings_bins_r = 1 
    IP_outer_rings_bins_phi = 20 

    # --------------------------------------------------------------------------------------
    # Outer Pixes sections - OP
    OP_r_min = IP_r_max
    OP_r_max = 345.
    OP_b_z_max = 380.
    OP_incl_z_max = 1100.
    OP_ring0_r_max = 205.
    OP_ring1_r_max = 255.

    # --------------------------------------------------------------------------------------
    # Pixels - P
    P_r_min = IP_r_min
    P_r_max = OP_r_max

    # --------------------------------------------------------------------------------------
    S_r_min = P_r_max
    S_r_max = ITk_central_r_max
    S_z_mid = 1400.
    S_ec_r_max = 990.


    # Augmenting the GeoModel sqlite
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES 
            (37, 
            'root', 
            'ITk', 
            'cyl;0.,{ITk_r_max},-{ITk_z_max},{ITk_z_max}',
            'children;NegSector,Central,PosSector',
            'z', 
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (36, 
            'leaf', 
            'ITk/NegSector', 
            'cyl;i,i,-{ITk_z_max},-{ITK_central_z_max}',
            '',
            '',
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (35, 
            'container', 
            'ITk/Central', 
            'cyl;i,i,-{ITK_central_z_max},{ITK_central_z_max}',
            'children;BeamPipe,Detectors,Outer',
            'r', 
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (33, 
            'leaf', 
            'ITk/PosSector',
            'cyl;,i,i,{ITK_central_z_max},{ITk_z_max}',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (0, 
            'leaf',
            'ITk/Central/BeamPipe',
            'cyl;i,{BP_r_max},i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (32, 
            'container',
            'ITk/Central/Detectors',
            'cyl;{IP_r_min},{ITk_central_r_max},i,i',
            'children;Pixels,Strips',
            'r', 
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (34, 
            'leaf', 
            'ITk/Central/Outer', 
            'cyl;{ITk_central_r_max},i,i,i',
            '',
            '',
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (30, 
            'container',
            'ITk/Central/Detectors/Pixels',
            'cyl;i,{P_r_max},i,i',
            'children;InnerPixels,OuterPixels',
            'r',
            '')
    """
    )
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES
            (31, 
            'container',
            'ITk/Central/Detectors/Strips',
            'cyl;{S_r_min},i,i,i',
            'children;NegSector,Barrel,PosSector',
            'z', 
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (20, 
            'container',
            'ITk/Central/Detectors/Pixels/InnerPixels',
            'cyl;i,{IP_r_max},i,i',
            'children;NegEndcap,Barrel,PosEndcap',
            'z',
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (25, 
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels',
            'cyl;{OP_r_min},i,i,i',
            'children;NegEndcap,NegInclined,Barrel,PosInclined,PosEndcap',
            'z',
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (27,
            'container',
            'ITk/Central/Detectors/Strips/NegSector',
            'cyl;i,i,i,-{S_z_mid}',
            'children;NegEndcap,OuterGap',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (18,
            'leaf',
            'ITk/Central/Detectors/Strips/Barrel',
            'cyl;i,i,-{S_z_mid},{S_z_mid}',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (29,
            'container',
            'ITk/Central/Detectors/Strips/PosSector',
            'cyl;i,i,{S_z_mid},i',
            'children;PosEndcap,OuterGap',
            'r', '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (17,
            'leaf',
            'ITk/Central/Detectors/Strips/NegSector/NegEndcap',
            'cyl;i,{S_ec_r_max},i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (26,
            'leaf',
            'ITk/Central/Detectors/Strips/NegSector/OuterGap',
            'cyl;{S_ec_r_max},i,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (19,
            'leaf',
            'ITk/Central/Detectors/Strips/PosSector/PosEndcap', 
            'cyl;i,{S_ec_r_max},i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (28,
            'leaf',
            'ITk/Central/Detectors/Strips/PosSector/OuterGap',
            'cyl;{S_ec_r_max},i,i,i',
            '',
            '', '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (21, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap', 
            'cyl;i,i,i,-{OP_incl_z_max}',
            'children;Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (22, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined',
            'cyl;i,i,-{OP_incl_z_max},-{OP_b_z_max}',
            'children;Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (10, 
            'leaf', 
            'ITk/Central/Detectors/Pixels/OuterPixels/Barrel',
            'cyl;i,i,-{OP_b_z_max},{OP_b_z_max}',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (23, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined',
            'cyl;i,i,{OP_b_z_max},{OP_incl_z_max}',
            'children;Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (24, 
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap',
            'cyl;i,i,{OP_incl_z_max},i',
            'children;Ring0,Ring1,Ring2',
            'r',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (4,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring0', 
            'cyl;i,{OP_ring0_r_max},i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (5,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring1',
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (6, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring2',
            'cyl;{OP_ring1_r_max},i,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES
            (7, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring0',
            'cyl;i,{OP_ring0_r_max},i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (8,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring1',
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (9,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring2',
            'cyl;{OP_ring1_r_max},i,i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (11,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring0',
            'cyl;i,{OP_ring0_r_max},i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (12,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring1',
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},i,i',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (13, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring2',
            'cyl;{OP_ring1_r_max},i,i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (14, 
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring0',
            'cyl;i,{OP_ring0_r_max},i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES
            (15,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring1',
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},i,i',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (16,
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring2',
            'cyl;{OP_ring1_r_max},i,i,i',
            '',
            '',
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (1, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/NegEndcap',
            'cyl;i,i,i,-{IP_b_z_max}',
            '',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (2,
            'container',
            'ITk/Central/Detectors/Pixels/InnerPixels/Barrel',
            'cyl;i,i,-{OP_ring1_r_max},{IP_b_z_max}',
            'children;Layer0,Layer1,Layer2,Layer3',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (200,
            'container',
            'ITk/Central/Detectors/Pixels/InnerPixels/Barrel',
            'cyl;i,i,-{OP_ring1_r_max},{IP_b_z_max}',
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
        f"""
    INSERT INTO Blueprint VALUES
            (3, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/PosEndcap',
            'cyl;i,i,{OP_ring1_r_max},i',
            '',
            '',
            '')
    """
    )

    connection.commit()

    print(">> Commited and done.")
