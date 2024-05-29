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
    print('>> Performing', cmd)
    os.system(cmd)  
    
    
    # Creating the database connection 
    connection = sqlite3.connect(args.output)
    cursor = connection.cursor()
    cursor.execute("CREATE TABLE Blueprint(id INT, type TEXT, name TEXT, bounds TEXT, internals TEXT, binnings TEXT, materials TEXT)")
    
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

    # Augmenting the GeoModel sqlite 
    cursor.execute("""
    INSERT INTO Blueprint VALUES
            (0, 'container', 'ITk', 'cyl;0.,1180,-3500.,3500.','ITk/NegEndplate,ITk/Central,ITk/PosEndplate','z', '')
    """)
    cursor.execute("""
    INSERT INTO Blueprint VALUES
            (100, 'leaf', 'ITk/NegEndplate', 'cyl;i,i,-3500.,-3050.','','', '')
    """)
    cursor.execute("""
    INSERT INTO Blueprint VALUES
            (101, 'leaf', 'ITk/Central', 'cyl;i,i,-3050.,3050.','','', '')
    """)
    cursor.execute("""
    INSERT INTO Blueprint VALUES
            (102, 'leaf', 'ITk/PosEndplate', 'cyl;,i,i,3050.,3500','','', '')
    """)
    
    connection.commit()
