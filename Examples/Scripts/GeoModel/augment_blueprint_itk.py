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

    volCount = 0

    volCount += 1
    cursor.execute(
        "CREATE TABLE Blueprint(id INT, type TEXT, name TEXT, bounds TEXT, internals TEXT, binnings TEXT, materials TEXT)"
    )

    # Bounds nomenclature is:
    # (type;[values])
    # type: cyl, box
    # values:
    #  - cyl, [rmin, rmax, zmin, zmax]
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
   
    IP_iec_z_max = 1090.

    # IP_iec z-positions: 265, 293, 324, 359, 398, 439, 488, 545, 606, 677, 751, 837, 927, 1028, 
    IP_iec_searchsplits = [ -IP_iec_z_max, -1000, -900, -800, -700, -620, -580, -520, -450, -410, -380, -340, -305, -280, -IP_b_z_max ]

    # IP_ec z-positions: 1105, 1144, 1231, 1274, 1361, 1405, 1505, 1555, 1667, 1723, 1848, 1911, 2122, 2359, 2623
    IP_ec_searchsplits =  [ -2700, -2500, -2200, -2000, -1880, -1800, -1700, -1600, -1520, -1450, -1380, -1300, -1250, -1200, -1120, -IP_iec_z_max ]

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
    OP_incl_z_max = IP_iec_z_max
    OP_ring0_r_max = 205.
    OP_ring1_r_max = 265.


    OP_b_searchsplits = [ OP_r_min, 220 , 280, OP_r_max ]

    OP_incl_ring0_searchsplits = [ -1090, -900, -750, -650, -550, -450, -375 ]
    OP_incl_ring1_searchsplits = [ -1090, -950, -850, -750, -650, -570, -500, -420, -375]
    OP_incl_ring2_searchsplits = [  -1090, -990, -870,-770,-700,-625, -545, -485, -420, -375] 

    OP_ec_ring0_searchsplits = [ -3000, -2800, -2500, -2350, -2100, -1900, -1700, -1550, -1400, -1300, -1200, -1100 ]
    OP_ec_ring1_searchsplits = [ -3000, -2600, -2300, -2000, -1800, -1600, -1400, -1200, -1100 ]
    OP_ec_ring2_searchsplits = [  -3000, -2700, -2400, -2100, -1900, -1700, -1500, -1300, -1200, -1100 ] 

    # Binnings
    # opBarrelBinnings = [ [ 18, 32], [18, 44], [18, 56] ]
    OP_b_z_bins = 18
    OP_b_phi_bins = [ 32, 44, 56 ]
    OP_incl_phi_bins = [ 32, 44, 56 ]
    OP_ec_phi_bins = [ 32, 44, 56 ]

    # --------------------------------------------------------------------------------------
    # Pixels - P
    P_r_min = IP_r_min
    P_r_max = OP_r_max

    # STRIPS --------------------------------------------------------------------------------
    S_r_min = P_r_max
    S_r_max = ITk_central_r_max
    S_z_mid = 1400.
    S_ec_r_max = 990.

    S_b_searchsplits = [ 350 , 500, 700, 850 , 1100] 
    S_ec_searchsplits = [ -3000 , -2700, -2400, -2100 , -1800, -1600, -1450] 


    # Strip binning: sBarrelBinnings = [ [ 56, 28*4], [56, 40*4], [56, 56*2], [56, 72*2] ]
    S_b_z_bins = 56
    S_b_phi_bins = [ 28, 40, 56, 72 ]


    # --------------------------------------------------------------------------------------

    volCount = 0

    volCount += 1    
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:110',
            'ITk',
            'cyl,0,{ITk_r_max},-{ITk_z_max},{ITk_z_max}',
            'children:BeamPipe,Container',
            'r', 
            '')
    """
    )

    volCount += 1    
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/BeamPipe',
            'cyl,e,{BP_r_max},e,e',
            '',
            '', 
            'p2:z,bound,1000,0,*,*;phi,closed,1,0,*,*')
    """
    )

    # Augmenting the GeoModel sqlite
    volCount += 1
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES 
            ({volCount}, 
            'root', 
            'ITk/Container', 
            'cyl,{BP_r_max},e,e,e',
            'children:NegSector,Central,PosSector',
            'z', 
            '')
    """
    )
    
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf', 
            'ITk/Container/NegSector', 
            'cyl,e,e,-{ITk_z_max},-{ITK_central_z_max}',
            '',
            '',
            '')
    """
    )
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:100', 
            'ITk/Container/Central', 
            'cyl,e,e,-{ITK_central_z_max},{ITK_central_z_max}',
            'children:Detectors,Outer',
            'r', 
            '')
    """
    )
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf', 
            'ITk/Container/PosSector',
            'cyl,,e,e,{ITK_central_z_max},{ITk_z_max}',
            '',
            '', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:98',
            'ITk/Container/Central/Detectors',
            'cyl,{IP_r_min},{ITk_central_r_max},e,e',
            'children:Pixels,Strips',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf', 
            'ITk/Container/Central/Outer', 
            'cyl,{ITk_central_r_max},e,e,e',
            '',
            '',
            '')
    """
    )
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:45',
            'ITk/Container/Central/Detectors/Pixels',
            'cyl,e,{P_r_max},e,e',
            'children:InnerPixels,OuterPixels',
            'r',
            '')
    """
    )
    volCount += 1
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:55',
            'ITk/Container/Central/Detectors/Strips',
            'cyl,{S_r_min},e,e,e',
            'children:NegSector,Barrel,PosSector',
            'z', 
            '')
    """
    )
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:15',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels',
            'cyl,e,{IP_r_max},e,e',
            'children:NegOuterEndcap,NegInnerEndcap,Barrel,PosInnerEndcap,PosOuterEndcap',
            'z',
            '')
    """
    )
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:14',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels',
            'cyl,{OP_r_min},e,e,e',
            'children:NegEndcap,NegInclined,Barrel,PosInclined,PosEndcap',
            'z',
            '')
    """
    )
    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:57',
            'ITk/Container/Central/Detectors/Strips/NegSector',
            'cyl,e,e,e,-{S_z_mid}',
            'children:NegEndcap,OuterGap',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:50',
            'ITk/Container/Central/Detectors/Strips/Barrel',
            'cyl,e,e,-{S_z_mid},{S_z_mid}',
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
        volCount += 1
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Strips/Barrel/Layer{i-1}',
            'cyl,i+2,i+2,e,e',
            'layer:kdt,cyl,{S_b_searchsplits[i-1]},{sp},e,e',
            'z,bound,{S_b_z_bins},1;phi,closed,{S_b_phi_bins[i-1]},1', 
            '')
    """
    )


    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:56',
            'ITk/Container/Central/Detectors/Strips/PosSector',
            'cyl,e,e,{S_z_mid},e',
            'children:PosEndcap,OuterGap',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:51',
            'ITk/Container/Central/Detectors/Strips/NegSector/NegEndcap',
            'cyl,e,{S_ec_r_max},e,e',
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
        volCount += 1
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Strips/NegSector/NegEndcap/Disk{len(S_ec_searchsplits)-1-i}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{S_ec_searchsplits[i-1]},{sp}',
            'r,bound,12,1;phi,closed,64,1', 
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Strips/NegSector/OuterGap',
            'cyl,{S_ec_r_max},e,e,e',
            '',
            '',
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:52',
            'ITk/Container/Central/Detectors/Strips/PosSector/PosEndcap', 
            'cyl,e,{S_ec_r_max},e,e',
            'children:*,Disk0,*,Disk1,*,Disk2,*,Disk3,*,Disk4,*,Disk5,*',
            'z',
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Strips/PosSector/OuterGap',
            'cyl,{S_ec_r_max},e,e,e',
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
        volCount += 1
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Strips/PosSector/PosEndcap/Disk{i-1}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{S_ec_searchsplits[i-1]},{sp}',
            'r,bound,12,1;phi,closed,64,1', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:81', 
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap', 
            'cyl,e,e,e,-{OP_incl_z_max}',
            'children:Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:82', 
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined',
            'cyl,e,e,-{OP_incl_z_max},-{OP_b_z_max}',
            'children:Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:83', 
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/Barrel',
            'cyl,e,e,-{OP_b_z_max},{OP_b_z_max}',
            'children:*,Layer0,*,Layer1,*,Layer2,*',
            'r',
            '')
    """
    )

    for i,sp in enumerate(OP_b_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        cursor.execute(
        f"""
        INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/Barrel/Layer{i-1}',
            'cyl,i+2,i+2,e,e',
            'layer:kdt,cyl,{OP_b_searchsplits[i-1]},{sp},e,e',
            'z,bound,{OP_b_z_bins},1;phi,closed,{OP_b_phi_bins[i-1]},1', 
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:84', 
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined',
            'cyl,e,e,{OP_b_z_max},{OP_incl_z_max}',
            'children:Ring0,Ring1,Ring2',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:85',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap',
            'cyl,e,e,{OP_incl_z_max},e',
            'children:Ring0,Ring1,Ring2',
            'r',
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:23',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring0', 
            'cyl,e,{OP_ring0_r_max},e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,ECRing8,*,ECRing9,*,ECRing10,*',
            'z',
            '')
    """
    )

    for i,sp in enumerate(OP_ec_ring0_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        cursor.execute(
        f"""INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring0/ECRing{len(OP_ec_ring0_searchsplits)-1-i}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{OP_ec_ring0_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_ec_phi_bins[0]},1', 
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:33',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring1',
            'cyl,{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,*',
            'z',
            '')
        """
        )

    for i,sp in enumerate(OP_ec_ring1_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        cursor.execute(
        f"""INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring1/ECRing{len(OP_ec_ring1_searchsplits)-1-i}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{OP_ec_ring1_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_ec_phi_bins[1]},1', 
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:43',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring2',
            'cyl,{OP_ring1_r_max},e,e,e',
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
        volCount += 1
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegEndcap/Ring2/ECRing{len(OP_ec_ring2_searchsplits)-1-i}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{OP_ec_ring2_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_ec_phi_bins[2]},1', 
            '')
        """
        )

    volCount += 1
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:21',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring0',
            'cyl,e,{OP_ring0_r_max},e,e',
            'children:*,InclRing0,*,InclRing1,*,InclRing2,*,InclRing3,*,InclRing4,*,InclRing5,*',
            'z', 
            '')
    """
    )

   # Inclined Layers
    for i,sp in enumerate(OP_incl_ring0_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != len(OP_incl_ring0_searchsplits)-1 else 'i+2,e'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring0/InclRing{len(OP_incl_ring0_searchsplits)-1-i}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{OP_incl_ring0_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_incl_phi_bins[0]},1',   
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:31',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring1',
            'cyl,{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'children:*,InclRing0,*,InclRing1,*,InclRing2,*,InclRing3,*,InclRing4,*,InclRing5,*,InclRing6,*,InclRing7,*',
            'z', 
            '')
    """
    )

   # Inclined Layers
    for i,sp in enumerate(OP_incl_ring1_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != len(OP_incl_ring1_searchsplits)-1 else 'i+2,e'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring1/InclRing{len(OP_incl_ring1_searchsplits)-1-i}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{OP_incl_ring1_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_incl_phi_bins[1]},1',   
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:41',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring2',
            'cyl,{OP_ring1_r_max},e,e,e',
            'children:*,InclRing0,*,InclRing1,*,InclRing2,*,InclRing3,*,InclRing4,*,InclRing5,*,InclRing6,*,InclRing7,*,InclRing8,*',
            'z',
            '')
    """
    )

   # Inclined Layers
    for i,sp in enumerate(OP_incl_ring2_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != len(OP_incl_ring2_searchsplits)-1 else 'i+2,e'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/NegInclined/Ring2/InclRing{len(OP_incl_ring2_searchsplits)-1-i}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{OP_incl_ring2_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_incl_phi_bins[2]},1',   
            '')
        """
        )


    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:24',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring0',
            'cyl,e,{OP_ring0_r_max},e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,ECRing8,*,ECRing9,*,ECRing10,*',
            'z',
            '')
    """
    )

    OP_ec_ring0_searchsplits = [ -1 * ss for ss in OP_ec_ring0_searchsplits ]
    OP_ec_ring0_searchsplits.reverse()
    for i,sp in enumerate(OP_ec_ring0_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        cursor.execute(
                f"""
        INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring0/ECRing{i-1}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{OP_ec_ring0_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_ec_phi_bins[0]},1',   
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:34',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring1',
            'cyl,{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'children:*,ECRing0,*,ECRing1,*,ECRing2,*,ECRing3,*,ECRing4,*,ECRing5,*,ECRing6,*,ECRing7,*',
            'z',
            '')
    """
    )
    OP_ec_ring1_searchsplits = [ -1 * ss for ss in OP_ec_ring1_searchsplits ]
    OP_ec_ring1_searchsplits.reverse()
    for i,sp in enumerate(OP_ec_ring1_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        cursor.execute(
                f"""
        INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring1/ECRing{i-1}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{OP_ec_ring1_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_ec_phi_bins[1]},1',   
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:44',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring2',
            'cyl,{OP_ring1_r_max},e,e,e',
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
        volCount += 1
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosEndcap/Ring2/ECRing{i-1}',
            'cyl,e,e,i+2,i+2',
            'layer:kdt,cyl,e,e,{OP_ec_ring2_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_ec_phi_bins[2]},1',   
            '')
        """
        )
        

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:22',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring0',
            'cyl,e,{OP_ring0_r_max},e,e',
            'children:*,InclRing0,*,InclRing1,*,InclRing2,*,InclRing3,*,InclRing4,*,InclRing5,*',
            'z', 
            '')
    """
    )

    OP_incl_ring0_searchsplits = [ -1 * ss for ss in OP_incl_ring0_searchsplits ]
    OP_incl_ring0_searchsplits.reverse()
    for i,sp in enumerate(OP_incl_ring0_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != 1 else 'e,i+2'
        cursor.execute(
        f"""           
      INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring0/InclRing{i-1}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{OP_incl_ring0_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_incl_phi_bins[0]},1',   
            '')
        """
        )

    volCount += 1
    cursor.execute(
       f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:32',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring1',
            'cyl,{OP_ring0_r_max},{OP_ring1_r_max},e,e',
            'children:*,InclRing0,*,InclRing1,*,InclRing2,*,InclRing3,*,InclRing4,*,InclRing5,*,InclRing6,*,InclRing7,*',
            'z', 
            '')
    """
    )

    # Inclined Layers
    OP_incl_ring1_searchsplits = [ -1 * ss for ss in OP_incl_ring1_searchsplits ]
    OP_incl_ring1_searchsplits.reverse()
    for i,sp in enumerate(OP_incl_ring1_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != 1 else 'e,i+2'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring1/InclRing{i-1}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{OP_incl_ring1_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_incl_phi_bins[1]},1',   
            '')
        """
        )


    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:42',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring2',
            'cyl,{OP_ring1_r_max},e,e,e',
            'children:*,InclRing0,*,InclRing1,*,InclRing2,*,InclRing3,*,InclRing4,*,InclRing5,*,InclRing6,*,InclRing7,*,InclRing8,*',
            'z',
            '')
    """
    )

    # Inclined Layers
    OP_incl_ring2_searchsplits = [ -1 * ss for ss in OP_incl_ring2_searchsplits ]
    OP_incl_ring2_searchsplits.reverse()
    
    for i,sp in enumerate(OP_incl_ring2_searchsplits):
           if i == 0:
                  continue
           # now make the layers
           volCount += 1
           zbest = 'i+2,i+2' if i != 1 else 'e,i+2'
           cursor.execute(
                f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/OuterPixels/PosInclined/Ring2/InclRing{i-1}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{OP_incl_ring2_searchsplits[i-1]},{sp}',
            'phi,closed,{OP_incl_phi_bins[2]},1',   
            '')
        """
        )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:13',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/NegOuterEndcap',
            'cyl,e,e,e,-{IP_iec_z_max}',
            'children:*,Disk0,*,Disk1,*,Disk2,*,Disk3,*,Disk4,*,Disk5,*,Disk6,*,Disk7,*,Disk8,*,Disk9,*,Disk10,*,Disk11,*,Disk12,*,Disk13,*',
            'z', 
            '')
    """
    )

    # Endcap Layers
    for i,sp in enumerate(IP_ec_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != len(IP_ec_searchsplits)-1 else 'i+2,e'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/NegOuterEndcap/Disk{len(IP_ec_searchsplits)-1-i}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{IP_ec_searchsplits[i-1]},{sp}',
            'phi,closed,30,1', 
            '')
        """
        )


    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:11',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/NegInnerEndcap',
            'cyl,e,e,-{IP_iec_z_max},-{IP_b_z_max}',
            'children:Disk0,*,Disk1,*,Disk2,*,Disk3,*,Disk4,*,Disk5,*,Disk6,*,Disk7,*,Disk8,*,Disk9,*,Disk10,*,Disk11,*,Disk12,*,Disk13,*',
            'z', 
            '')
    """
    )

    # Endcap Disks
    for i,sp in enumerate(IP_iec_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != len(IP_iec_searchsplits)-1 else 'i+2,e'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/NegInnerEndcap/Disk{len(IP_iec_searchsplits)-1-i}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{IP_iec_searchsplits[i-1]},{sp}',
            'r,bound,2,1;phi,closed,20,1', 
            '')
        """
        )  


    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:10',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/Barrel',
            'cyl,e,e,-{IP_b_z_max},{IP_b_z_max}',
            'children:Layer0,*,Layer1,*',
            'r', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/Barrel/Layer0',
            'cyl,i+2,i+2,e,e',
            'layer:kdt,cyl,30,80,e,e',
            'z,bound,12,1;phi,closed,12,1', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/Barrel/Layer1',
            'cyl,i+2,i+2,e,e',
            'layer:kdt,cyl,80,120,e,e',
            'z,bound,12,1;phi,closed,20,1', 
            '')
    """
    )

    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:12',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/PosInnerEndcap',
            'cyl,e,e,{IP_b_z_max},{IP_iec_z_max}',
           'children:Disk0,*,Disk1,*,Disk2,*,Disk3,*,Disk4,*,Disk5,*,Disk6,*,Disk7,*,Disk8,*,Disk9,*,Disk10,*,Disk11,*,Disk12,*,Disk13,*',
            'z', 
            '')
    """
    )

    # Endcap Disks
    IP_iec_searchsplits = [ -1 * ss for ss in IP_iec_searchsplits ]
    IP_iec_searchsplits.reverse()
    for i,sp in enumerate(IP_iec_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != 1 else 'e,i+2'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/PosInnerEndcap/Disk{i-1}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{IP_iec_searchsplits[i-1]},{sp}',
            'r,bound,2,1;phi,closed,20,1',  
            '')
        """
        )



    volCount += 1
    cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'container:14',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/PosOuterEndcap',
            'cyl,e,e,{IP_iec_z_max},e',
             'children:*,Disk0,*,Disk1,*,Disk2,*,Disk3,*,Disk4,*,Disk5,*,Disk6,*,Disk7,*,Disk8,*,Disk9,*,Disk10,*,Disk11,*,Disk12,*,Disk13,*',
            'z', 
            '')
    """
    )

    # Endcap layers
    IP_ec_searchsplits = [ -1 * ss for ss in IP_ec_searchsplits ]
    IP_ec_searchsplits.reverse()
    for i,sp in enumerate(IP_ec_searchsplits):
        if i == 0:
            continue
        # now make the layers
        volCount += 1
        zbest = 'i+2,i+2' if i != 1 else 'e,i+2'
        cursor.execute(
        f"""
    INSERT INTO Blueprint VALUES
            ({volCount}, 
            'leaf',
            'ITk/Container/Central/Detectors/Pixels/InnerPixels/PosOuterEndcap/Disk{i-1}',
            'cyl,e,e,{zbest}',
            'layer:kdt,cyl,e,e,{IP_ec_searchsplits[i-1]},{sp}',
            'phi,closed,30,1', 
            '')
        """
        )


    connection.commit()

    print(">> Commited and done.")
