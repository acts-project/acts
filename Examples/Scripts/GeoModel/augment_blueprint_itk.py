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
    # type: cyl| box
    # values:
    #  - cyl| [rmin, rmax, zmin, zmax]
    #  - box: [xmin, xmax, ymin, ymax, zmin, zmax]
    # special characters: e ... external (from mother),
    #                     c ... calculated (from container siblings),
    #                     i ... internal (from internals) with default value
    #                     iX ... internal (from internals) with value X
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
    IP_r_max = 130.
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
    OP_b_z_max = 376.
    OP_incl_z_max = 1100.
    OP_ring0_r_max = 205.
    OP_ring1_r_max = 265.

    OP_b_searchsplits = [ OP_r_min, 220 , 280, OP_r_max ]

    OP_incl_ring2_searchsplits = [ ]

    OP_ec_ring0_searchsplits = [ ]
    OP_ec_ring1_searchsplits = [ -3000, -2600, -2300, -2000, -1800, -1600, -1400, -1200, -1100 ]
    OP_ec_ring2_searchsplits = [  -3000, -2700, -2400, -2100, -1900, -1700, -1500, -1300, -1200, -1100 ] 

    # --------------------------------------------------------------------------------------
    # Pixels - P
    P_r_min = IP_r_min
    P_r_max = OP_r_max

    # --------------------------------------------------------------------------------------
    S_r_min = P_r_max
    S_r_max = ITk_central_r_max
    S_z_mid = 1400.
    S_ec_r_max = 990.

    S_b_searchsplits = [ 350 , 500, 700, 850 , 1100] 
    S_ec_searchsplits = [ -3000 , -2700, -2400, -2100 , -1800, -1600, -1450] 


    # Augmenting the GeoModel sqlite
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES 
            (37, 
            'root', 
            'ITk', 
            'cyl|0.,{ITk_r_max},-{ITk_z_max},{ITk_z_max}',
            'children:NegSector,Central,PosSector',
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
            'cyl|e,e,-{ITk_z_max},-{ITK_central_z_max}',
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
            'cyl|e,e,-{ITK_central_z_max},{ITK_central_z_max}',
            'children:BeamPipe,Detectors,Outer',
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
            'cyl|,e,e,{ITK_central_z_max},{ITk_z_max}',
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
            'cyl|e,{BP_r_max},e,e',
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
            'cyl|{IP_r_min},{ITk_central_r_max},e,e',
            'children:Pixels,Strips',
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
            'cyl|{ITk_central_r_max},e,e,e',
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
            'cyl|e,{P_r_max},e,e',
            'children:InnerPixels,OuterPixels',
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
            'cyl|{S_r_min},e,e,e',
            'children:NegSector,Barrel,PosSector',
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
            'cyl|e,{IP_r_max},e,e',
            'children:NegOuterEndcap,NegInnerEndcap,Barrel,PosInnerEndcap,PosOuterEndcap',
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
            'cyl|{OP_r_min},e,e,e',
            'children:NegEndcap,NegInclined,Barrel,PosInclined,PosEndcap',
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
            'cyl|e,e,e,-{S_z_mid}',
            'children:NegEndcap,OuterGap',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (18,
            'container',
            'ITk/Central/Detectors/Strips/Barrel',
            'cyl|e,e,-{S_z_mid},{S_z_mid}',
            'children:*,Layer0,*,Layer1,*,Layer2,*,Layer3,*',
            'r',
            '')
    """
    )

    # Barrel Layers
    for i,sp in enumerate(S_b_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (2000{i},
            'leaf',
            'ITk/Central/Detectors/Strips/Barrel/Layer{i-1}',
            'cyl|i+2,i+2,e,e',
            'layer:kdt|cyl|{S_b_searchsplits[i-1]},{sp},e,e',
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
            'cyl|e,e,{S_z_mid},e',
            'children:PosEndcap,OuterGap',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (17,
            'container',
            'ITk/Central/Detectors/Strips/NegSector/NegEndcap',
            'cyl|e,{S_ec_r_max},e,e',
            'children:*,Disk5,*,Disk4,*,Disk3,*,Disk2,*,Disk1,*,Disk0,*',
            'z',
            '')
    """
    )

   # Endcap Layers
    for i,sp in enumerate(S_ec_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (3000{i},
            'leaf',
            'ITk/Central/Detectors/Strips/NegSector/NegEndcap/Disk{len(S_ec_searchsplits)-1-i}',
            'cyl|e,e,i+2,i+2',
            'layer:kdt|cyl|e,e,{S_ec_searchsplits[i-1]},{sp}',
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
            'cyl|{S_ec_r_max},e,e,e',
            '',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (19,
            'container',
            'ITk/Central/Detectors/Strips/PosSector/PosEndcap', 
            'cyl|e,{S_ec_r_max},e,e',
            'children:*,Disk0,*,Disk1,*,Disk2,*,Disk3,*,Disk4,*,Disk5,*',
            'z',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (28,
            'leaf',
            'ITk/Central/Detectors/Strips/PosSector/OuterGap',
            'cyl|{S_ec_r_max},e,e,e',
            '',
            '', 
            '')
    """
    )

   # Barrel Layers
    S_ec_searchsplits = [ -1 * ss for ss in S_ec_searchsplits ]
    S_ec_searchsplits.reverse()
    for i,sp in enumerate(S_ec_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (4000{i},
            'leaf',
            'ITk/Central/Detectors/Strips/PosSector/PosEndcap/Disk{i-1}',
            'cyl|e,e,i+2,i+2',
            'layer:kdt|cyl|e,e,{S_ec_searchsplits[i-1]},{sp}',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (21, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap', 
            'cyl|e,e,e,-{OP_incl_z_max}',
            'children:Ring0,Ring1,Ring2',
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
            'cyl|e,e,-{OP_incl_z_max},-{OP_b_z_max}',
            'children:Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (10, 
            'container', 
            'ITk/Central/Detectors/Pixels/OuterPixels/Barrel',
            'cyl|e,e,-{OP_b_z_max},{OP_b_z_max}',
            'children:*,Layer0,*,Layer1,*,Layer2,*',
            'r',
            '')
    """
    )

    for i,sp in enumerate(OP_b_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""
        INSERT INTO Blueprint VALUES
            (1000{i+1},
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/Barrel/Layer{i-1}',
            'cyl|i+2,i+2,e,e',
            'layer:kdt|cyl|{OP_b_searchsplits[i-1]},{sp},e,e',
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
            'cyl|e,e,{OP_b_z_max},{OP_incl_z_max}',
            'children:Ring0,Ring1,Ring2',
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
            'cyl|e,e,{OP_incl_z_max},e',
            'children:Ring0,Ring1,Ring2',
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
            'cyl|e,{OP_ring0_r_max},e,e',
            'layer:kdt|cyl|e,{OP_ring0_r_max},e,e',
            '',
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (5,
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring1',
            'cyl|{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,*',
            'z',
            '')
    """
    )

    for i,sp in enumerate(OP_ec_ring1_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""INSERT INTO Blueprint VALUES
            (5008{i},
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring1/ECRing{len(OP_ec_ring1_searchsplits)-1-i}',
            'cyl|e,e,i+2,i+2',
            'layer:kdt|cyl|e,e,{OP_ec_ring1_searchsplits[i-1]},{sp}',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (6, 
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring2',
            'cyl|{OP_ring1_r_max},e,e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,*,ECRing8,*',
            'z', 
            '')
    """
    )

    # Endcap Layers
    for i,sp in enumerate(OP_ec_ring2_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (5009{i},
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring2/ECRing{len(OP_ec_ring2_searchsplits)-1-i}',
            'cyl|e,e,i+2,i+2',
            'layer:kdt|cyl|e,e,{OP_ec_ring2_searchsplits[i-1]},{sp}',
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
            'cyl|e,{OP_ring0_r_max},e,e',
            'layer:kdt|cyl|e,{OP_ring0_r_max},e,e',
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
            'cyl|{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'layer:kdt|cyl|{OP_ring0_r_max},{OP_ring1_r_max},e,e',
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
            'cyl|{OP_ring1_r_max},e,e,e',
            'layer:kdt|cyl|{OP_ring1_r_max},e,e,e',
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
            'cyl|e,{OP_ring0_r_max},e,e',
            'layer:kdt|cyl|e,{OP_ring0_r_max},e,e',
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
            'cyl|{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,*',
            '',
            '')
    """
    )
    OP_ec_ring1_searchsplits = [ -1 * ss for ss in OP_ec_ring1_searchsplits ]
    OP_ec_ring1_searchsplits.reverse()
    for i,sp in enumerate(OP_ec_ring1_searchsplits):
                if i == 0:
                        continue
                # now make the layers
                cursor.execute(
                f"""
INSERT INTO Blueprint VALUES
            (6000{i},
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring1/ECRing{i-1}',
            'cyl|e,e,i+2,i+2',
            'layer:kdt|cyl|e,e,{OP_ec_ring1_searchsplits[i-1]},{sp}',
            '', 
            '')
    """
    )

    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (13, 
            'container',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring2',
            'cyl|{OP_ring1_r_max},e,e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,*,ECRing8,*',
            'z', 
            '')
    """
    )

    # Endcap Layers
    OP_ec_ring2_searchsplits = [ -1 * ss for ss in OP_ec_ring2_searchsplits ]
    OP_ec_ring2_searchsplits.reverse()
    for i,sp in enumerate(OP_ec_ring2_searchsplits):
        if i == 0:
            continue
        # now make the layers
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            (6000{i},
            'leaf',
            'ITk/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring2/ECRing{i-1}',
            'cyl|e,e,i+2,i+2',
            'layer:kdt|cyl|e,e,{OP_ec_ring2_searchsplits[i-1]},{sp}',
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
            'cyl|e,{OP_ring0_r_max},e,e',
            'layer:kdt|cyl|e,{OP_ring0_r_max},e,e',
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
            'cyl|{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'layer:kdt|cyl|{OP_ring0_r_max},{OP_ring1_r_max},e,e',
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
            'cyl|{OP_ring1_r_max},e,e,e',
            'layer:kdt|cyl|{OP_ring1_r_max},e,e,e',
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
            'cyl|e,e,e,-{IP_iec_z_max}',
            'layer:kdt|cyl|e,e,e,-{IP_iec_z_max}',
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
            'cyl|e,e,-{IP_iec_z_max},-{IP_b_z_max}',
            'layer:kdt|cyl|e,e,-{IP_iec_z_max},-{IP_b_z_max}',
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
            'cyl|e,e,-{IP_b_z_max},{IP_b_z_max}',
            'children:Layer0,*,Layer1,*',
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
            'cyl|i+2,i+2,e,e',
            'layer:kdt|cyl|30,80,e,e',
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
            'cyl|i+2,i+2,e,e',
            'layer:kdt|cyl|80,120,e,e',
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
            'cyl|e,e,{IP_b_z_max},{IP_iec_z_max}',
            'layer:kdt|cyl|e,e,{IP_b_z_max},{IP_iec_z_max}',
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
            'cyl|e,e,{IP_iec_z_max},e',
            'layer:kdt|cyl|e,e,{IP_iec_z_max},e',
            '',
            '')
    """
    )

    connection.commit()

    print(">> Commited and done.")
