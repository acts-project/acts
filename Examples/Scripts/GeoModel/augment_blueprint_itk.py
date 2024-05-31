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
    IP_b_z_max = 255.
   
    IP_iec_z_max = 1100.

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
            'cyl;e,e,-{ITk_z_max},-{ITK_central_z_max}',
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
            'cyl;e,e,-{ITK_central_z_max},{ITK_central_z_max}',
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
            'cyl;,e,e,{ITK_central_z_max},{ITk_z_max}',
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
            'cyl;e,{BP_r_max},e,e',
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
            'cyl;{IP_r_min},{ITk_central_r_max},e,e',
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
            'cyl;{ITk_central_r_max},e,e,e',
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
            'cyl;e,{P_r_max},e,e',
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
            'cyl;{S_r_min},e,e,e',
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
            'cyl;e,{IP_r_max},e,e',
            'children;NegOuterEndcap,NegInnerEndcap,Barrel,PosInnerEndcap,PosOuterEndcap',
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
            'cyl;{OP_r_min},e,e,e',
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
            'cyl;e,e,e,-{S_z_mid}',
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
            'cyl;e,e,-{S_z_mid},{S_z_mid}',
            'layer;kdt;cyl;e,e,-{S_z_mid},{S_z_mid}',
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
            'cyl;e,e,{S_z_mid},e',
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
            'cyl;e,{S_ec_r_max},e,e',
            'layer;kdt;cyl;e,{S_ec_r_max},e,e',
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
            'cyl;{S_ec_r_max},e,e,e',
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
            'cyl;e,{S_ec_r_max},e,e',
            'layer;kdt;cyl;e,{S_ec_r_max},e,e',
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
            'cyl;{S_ec_r_max},e,e,e',
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
            'cyl;e,e,e,-{OP_incl_z_max}',
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
            'cyl;e,e,-{OP_incl_z_max},-{OP_b_z_max}',
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
            'cyl;e,e,-{OP_b_z_max},{OP_b_z_max}',
            'layer;kdt;cyl;e,e,-{OP_b_z_max},{OP_b_z_max}',
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
            'cyl;e,e,{OP_b_z_max},{OP_incl_z_max}',
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
            'cyl;e,e,{OP_incl_z_max},e',
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
            'cyl;e,{OP_ring0_r_max},e,e',
            'layer;kdt;cyl;e,{OP_ring0_r_max},e,e',
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
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'layer;kdt;cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
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
            'cyl;{OP_ring1_r_max},e,e,e',
            'layer;kdt;cyl;{OP_ring1_r_max},e,e,e',
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
            'cyl;e,{OP_ring0_r_max},e,e',
            'layer;kdt;cyl;e,{OP_ring0_r_max},e,e',
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
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'layer;kdt;cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
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
            'cyl;{OP_ring1_r_max},e,e,e',
            'layer;kdt;cyl;{OP_ring1_r_max},e,e,e',
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
            'cyl;e,{OP_ring0_r_max},e,e',
            'layer;kdt;cyl;e,{OP_ring0_r_max},e,e',
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
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'layer;kdt;cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
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
            'cyl;{OP_ring1_r_max},e,e,e',
            'layer;kdt;cyl;{OP_ring1_r_max},e,e,e',
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
            'cyl;e,{OP_ring0_r_max},e,e',
            'layer;kdt;cyl;e,{OP_ring0_r_max},e,e',
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
            'cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'layer;kdt;cyl;{OP_ring0_r_max},{OP_ring1_r_max},e,e',
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
            'cyl;{OP_ring1_r_max},e,e,e',
            'layer;kdt;cyl;{OP_ring1_r_max},e,e,e',
            '',
            '')
    """
    )
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (1001, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/NegOuterEndcap',
            'cyl;e,e,e,-{IP_iec_z_max}',
            'layer;kdt;cyl;e,e,e,-{IP_iec_z_max}',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (1002, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/NegInnerEndcap',
            'cyl;e,e,-{IP_iec_z_max},-{IP_b_z_max}',
            'layer;kdt;cyl;e,e,-{IP_iec_z_max},-{IP_b_z_max}',
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
            'cyl;e,e,-{IP_b_z_max},{IP_b_z_max}',
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
            'cyl;e,i+2,e,e',
            'layer;kdt;cyl;e,50,e,e',
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
            'cyl;i+2,i+2,e,e',
            'layer;kdt;cyl;80,120,e,e',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (103, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/PosInnerEndcap',
            'cyl;e,e,{IP_b_z_max},{IP_iec_z_max}',
            'layer;kdt;cyl;e,e,{IP_b_z_max},{IP_iec_z_max}',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (104, 
            'leaf',
            'ITk/Central/Detectors/Pixels/InnerPixels/PosOuterEndcap',
            'cyl;e,e,{IP_iec_z_max},e',
            'layer;kdt;cyl;e,e,{IP_iec_z_max},e',
            '',
            '')
    """
    )

    connection.commit()

    print(">> Commited and done.")
