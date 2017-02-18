#!/usr/bin/env python
"""
Created on Jul 24, 2012
@author: maechlin
$Id: gen_site_corr_list.py 1624 2016-04-08 18:36:46Z fsilva $
"""

import sys
import sqlite3
from optparse import OptionParser
from load_eventfile import periods

def main():
    #print "Starting"
    sqlitedb = None
    eqname = None
    usage = "usage: %s [options]" % (sys.argv[0])
    parser = OptionParser(usage)
    parser.add_option("-d", "--db_file", dest="sqlitedb",
                      help="File name containing SQLite3 DB")    
    parser.add_option("-e", "--event_name", dest="eqname",
                      help="Specify eqname for desired stations")
    (options, args) = parser.parse_args()

    # 
    # Point to user supplied DB is needed
    #
    if options.sqlitedb is not None:
        sqlitedb = options.sqlitedb
    if options.eqname is not None:
        eqname = options.eqname

    if sqlitedb is None or eqname is None:
        print "Must provide both sqlitedb and event name!"
        sys.exit(1)
 
    #print "Starting query"
    # Begin DB Query
    # Create a connection to the database.
    conn = sqlite3.connect(sqlitedb)
    # Create a cursor object to do the interacting.
    cursor = conn.cursor()

    # Get eqid from event name
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

    rows2 = cursor.execute('Select stationID, slat, slong, PGAr, '
                           'p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, '
                           'p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, '
                           'p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, '
                           'p30, p31, p32, p33, p34, p35, p36, p37, p38, p39, '
                           'p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, '
                           'p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, '
                           'p60, p61, p62 from RSNs where eqid = ? '
                           'order by stationID asc', [eqid])
    # Print header
    print "#BBP Station Site Correction List for EQID=%d (%s)" % (eqid, eqname)
    print ('#StaName StaLat  StaLon  PGAr(g) '
           '%s %s %s %s %s %s %s %s %s %s '
           '%s %s %s %s %s %s %s %s %s %s '
           '%s %s %s %s %s %s %s %s %s %s '
           '%s %s %s %s %s %s %s %s %s %s '
           '%s %s %s %s %s %s %s %s %s %s '
           '%s %s %s %s %s %s %s %s %s %s '
           '%s %s %s' %
           (periods[0], periods[1], periods[2], periods[3],
            periods[4], periods[5], periods[6], periods[7],
            periods[8], periods[9], periods[10], periods[11],
            periods[12], periods[13], periods[14], periods[15],
            periods[16], periods[17], periods[18], periods[19],
            periods[20], periods[21], periods[22], periods[23],
            periods[24], periods[25], periods[26], periods[27],
            periods[28], periods[29], periods[30], periods[31],
            periods[32], periods[33], periods[34], periods[35],
            periods[36], periods[37], periods[38], periods[39],
            periods[40], periods[41], periods[42], periods[43],
            periods[44], periods[45], periods[46], periods[47],
            periods[48], periods[49], periods[50], periods[51],
            periods[52], periods[53], periods[54], periods[55],
            periods[56], periods[57], periods[58], periods[59],
            periods[60], periods[61], periods[62]))
    # Print results from database
    for row in rows2:
        print ('%s %7.3f %6.3f % 6.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f %06.4f '
               '%06.4f %06.4f %06.4f' %
               (row[0], row[1], row[2], row[3], row[4], row[5], row[6],
                row[7], row[8], row[9], row[10], row[11], row[12], row[13],
                row[14], row[15], row[16], row[17], row[18], row[19], row[20],
                row[21], row[22], row[23], row[24], row[25], row[26], row[27],
                row[28], row[29], row[30], row[31], row[32], row[33], row[34],
                row[35], row[36], row[37], row[38], row[39], row[40], row[41],
                row[42], row[43], row[44], row[45], row[46], row[47], row[48],
                row[49], row[50], row[51], row[52], row[53], row[54], row[55],
                row[56], row[57], row[58], row[59], row[60], row[61], row[62],
                row[63], row[64], row[65], row[66]))

    conn.close()

if __name__ == "__main__":
    main()
