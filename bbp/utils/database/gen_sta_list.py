#!/usr/bin/env python
"""
Created on Jul 24, 2012
@author: maechlin
$Id: gen_sta_list.py 1624 2016-04-08 18:36:46Z fsilva $
"""

import sys
import sqlite3
from optparse import OptionParser

def main():
    """
    Get station list from the database
    """

    sqlitedb = None
    eqname = None
    ucsb = False
    ucsb_vs30 = False
    usage = "usage: %s [options]" % (sys.argv[0])
    parser = OptionParser(usage)
    parser.add_option("-d", "--db_file", dest="sqlitedb",
                      help="File name containing SQLite3 DB")    
    parser.add_option("-e", "--event_name", dest="eqname",
                      help="Specify eqname for desired stations")
    parser.add_option("--ucsb", action="store_true", dest="ucsb",
                      help="Generate station list for the UCSB code")
    parser.add_option("--ucsb-vs30", action="store_true", dest="ucsb_vs30",
                      help="Generate vs30 station list for the UCSB code")
    (options, args) = parser.parse_args()

    # 
    # Point to user supplied DB is needed
    #
    if options.sqlitedb is not None:
        sqlitedb = options.sqlitedb
    if options.eqname is not None:
        eqname = options.eqname
    if options.ucsb is not None:
        ucsb = True
    if options.ucsb_vs30 is not None:
        ucsb_vs30 = True

    if sqlitedb is None or eqname is None:
        print "Must provide both sqlitedb and event name!"
        sys.exit(1)

    if ucsb and ucsb_vs30:
        print "Use only one ucsb option at a time!"
        sys.exit(1)

    # Begin DB Query
    # Create a connection to the database.
    conn = sqlite3.connect(sqlitedb)
    # Create a cursor object to do the interacting.
    cursor = conn.cursor()
    # Get eqid from event_name
    cursor.execute("select eqid from EVENTS where Event_Name = ?",
                   [eqname])
    eqid = cursor.fetchone()

    if eqid is None:
        print "Event %s not found in db" % (eqname)
        sys.exit(1)

    eqid = eqid[0]

    # Grab the columns with the time-zoned dates.
    #rows = c.execute('SELECT * FROM Events')
    #for row in rows:
    #    print row[0],row[1]

    # I think this query only works with one event in the db. if you
    # add a second, you will need to update the query to select by
    # eventid.
    #
    rows2 = cursor.execute("Select stationid, slat, slong, vs30, "
                           "Lowpass_Freq, Highpass_Freq, n_file "
                           "from RSNs where eqid = ?",
                           [eqid])
    if ucsb:
        # Generate UCSB station lists
        station_list = []
        # Create list of stations
        for row in rows2:
            station_list.append((row[2], row[1], row[0]))
        print len(station_list)
        for station in station_list:
            print ("%7.4f %7.4f %s" %
                   (station[0], station[1], station[2]))
    elif ucsb_vs30:
        # Generate the UCSB VS30 station list
        for row in rows2:
            print "%s %d" % (row[0], row[3])
    else:
        # Generate station list
        print "#BBP Station List for EQID=%d (%s)" % (eqid, eqname)
        print ("#SLong    SLat     RSN   Vs30(m/s) "
               "LoPass_Freq(Hz) HiPass_Freq(Hz)  Obs_File")
        for row in rows2:
            print ("%7.3f %6.3f  %s %5d %5.4f  %5.4f    %s" %
                   (row[2], row[1], row[0], row[3],
                    row[4], row[5], row[6]))
    conn.close()

if __name__ == "__main__":
    main()
