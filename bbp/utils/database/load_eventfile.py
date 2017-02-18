#!/usr/bin/env python
"""
Created on Jul 24, 2012
@author: maechlin
$Id: load_eventfile.py 1625 2016-04-08 18:41:45Z fsilva $
"""

from optparse import OptionParser
import sqlite3
import sys

periods = (0.010000,
           0.011000,
           0.012000,
           0.013000,
           0.015000,
           0.017000,
           0.020000,
           0.022000,
           0.025000,
           0.029000,
           0.032000,
           0.035000,
           0.040000,
           0.045000,
           0.050000,
           0.055000,
           0.060000,
           0.065000,
           0.075000,
           0.085000,
           0.100000,
           0.110000,
           0.120000,
           0.130000,
           0.150000,
           0.170000,
           0.200000,
           0.220000,
           0.240000,
           0.260000,
           0.280000,
           0.300000,
           0.350000,
           0.400000,
           0.450000,
           0.500000,
           0.550000,
           0.600000,
           0.650000,
           0.750000,
           0.850000,
           1.000000,
           1.100000,
           1.200000,
           1.300000,
           1.500000,
           1.700000,
           2.000000,
           2.200000,
           2.400000,
           2.600000,
           2.800000,
           3.000000,
           3.500000,
           4.000000,
           4.400000,
           5.000000,
           5.500000,
           6.000000,
           6.500000,
           7.500000,
           8.500000,
           10.000000)

def indx2period(indx):
    print "Index: %d is period: %f " % (periods[indx])

def main():
    sqlitedb = None
    eventfile = None
    eqname = None
    usage = "usage: %s [options]" % (sys.argv[0])
    parser = OptionParser(usage)
    
    parser.add_option("-d", "--db-file", dest="sqlitedb",
                      help="File containing SQLite3 DB")
    parser.add_option("-f", "--event-file", dest="eventfile",
                      help=("File containing event information "
                            "for the validation study"))
    parser.add_option("-n", "--eq-name", dest="eqname",
                      help="Event name (i.e. Northridge, Landers, etc")
    #
    # Start of parsing input parameters to db query
    #
    (options, args) = parser.parse_args()

    # 
    # Point to user supplied DB is needed
    #
    if options.sqlitedb is not None:
        sqlitedb = options.sqlitedb
        print "Querying BBP_DB file: %s" % (sqlitedb)
    if options.eventfile is not None:
        eventfile = options.eventfile
        print "Loading Event File: %s" % (eventfile)
    if options.eqname is not None:
        eqname = options.eqname
        print "EQID: %s" % (eqname)

    if sqlitedb is None or eventfile is None or eqname is None:
        print "Please provide event name and both database and event files!"
        sys.exit(1)

    lins = []
    f = open(eventfile, 'U') # This Universal file open fixes DOS versus Mac cr\lf issue
    for line in f:
        # print line
        res = line.split('\n')
        lins.append(res[0])

    lnum = 0
    for l in lins:
        if lnum == 0:  # First line in column header line
            lnum = lnum + 1
        elif lnum == 1: # Second line is first line of data. Get Event ID and check only this ID is in file
            lnum = lnum + 1
            fields = l.split()
            if len(fields) != 79: # Confirm expected number of columns
                print "Unexpected columns in file: %d" % (len(fields))
                sys.exit(1)

            #if int(fields[0]) != eqid: # Compare new eqids from first row to all other rows. Only one eqid expected in file
            #    #print "Multiple event IDS found in Event File"
            #    print "Mismatch Event IDs between Event ID: %d and Event ID: %d"%(eqid,int(fields[0]))
                
    print "Total Station Entries in Event File: %s is %d" % (eventfile, lnum-1)

    #
    # Confirm that captured eqid is already in event table
    #
    connection = sqlite3.connect(sqlitedb)
    cursor = connection.cursor()
    cursor.execute("select eqid from EVENTS where Event_Name = ?",
                   [eqname]) # Found this cursor format at url
    #http://stackoverflow.com/questions/228912/sqlite-parameter-substitution-problem for more details
    #
    eqid = cursor.fetchone()

    if eqid is None:
        print "Event %s not found in db" % (eqname)
        sys.exit(1)

    eqid = eqid[0]
    print "Event %s has eqid %d" % (eqname, eqid)

    #
    period_start_col = 16

    #
    # Write each station information to RSN table
    #
    lnum = 0
    for l in lins:
        if lnum == 0:  # First line in column header line
            lnum = lnum + 1
        else:
            lnum = lnum + 1
            fields = l.split()
            if len(fields) != 79: # Confirm expected number of columns
                print "Unexpected columns in file: %d" % (len(fields))
                exit(1)
            
            #eqid = 118 # Warning. Hard coded EQID for LOMAP event. 
            # it is set at the start of the file. Must change to pass as param
            recid = int(fields[0])
            staid = fields[1]
            study_name  = "BVS1" #TODO move this hard coded value to a config file or header file
            slat = float(fields[2]) #Stands for BBP Validation Study 1
            slong = float(fields[3])
            rrup = float(fields[4])
            vs30 = float(fields[5])
            hp_freq = 1.0/float(fields[6])
            lp_freq = 1.0/float(fields[7])
            n_file = fields[12]
            e_file = fields[13]
            z_file = fields[14]
            pgar = float(fields[15])
            #print eqid
            #print rsn
            #print n_file
            #print e_file
            #print z_file
            #print study_name
            #print slat
            #print slong
            #print rrup
            #print vs30
            #column order for rsn table
            #RSN
            #EQID
            #Study_Name
            #SLAT
            #SLONG
            #Rrup
            #VS30
            #n_file
            #e_file
            #z_file
            #cmd = 'insert into RSNs (RECID, EQID, StationID, Study_Name, SLAT, SLONG, Rrup,VS30, Lowpass_Freq, Highpass_Freq, n_file, e_file, z_file) values (%d,%d,"%s","%s",%f,%f,%f,%f,%f,%f,"%s","%s","%s")'%(recid,eqid,staid,study_name,slat,slong,rrup,vs30, hp_freq, lp_freq, n_file, e_file, z_file)
                                                         
            #print "Query str:%s"%(cmd)
            #cursor.execute(cmd)            
            #connection.commit()
    
            #
            #
            #
            #for indx,p in enumerate(periods):
            #    cmd = 'update RSNs (RECID,EQID,?'
            #
            cmd = ('insert into RSNs (RECID, EQID, StationID, Study_Name, '
                   'SLAT, SLONG, Rrup, VS30, Lowpass_Freq, '
                   'Highpass_Freq, n_file, e_file, z_file, PGAr, '
                   'p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, '
                   'p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, '
                   'p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, '
                   'p35, p36, p37, p38, p39, p40, p41, p42, p43, p44, p45, '
                   'p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, '
                   'p57, p58, p59, p60, p61, p62) values (%d, %d, "%s","%s", '
                   '%f, %f, %f, %f, %f, %f, "%s", "%s", "%s", %f, %f, %f, %f, '
                   '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, '
                   '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, '
                   '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, '
                   '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, '
                   '%f, %f, %f, %f)' %
                   (recid,
                    eqid, staid, study_name, slat, slong, rrup, vs30,
                    lp_freq, hp_freq, n_file, e_file, z_file, pgar,
                    float(fields[period_start_col]),
                    float(fields[period_start_col+1]),
                    float(fields[period_start_col+2]),
                    float(fields[period_start_col+3]),
                    float(fields[period_start_col+4]),
                    float(fields[period_start_col+5]),
                    float(fields[period_start_col+6]),
                    float(fields[period_start_col+7]),
                    float(fields[period_start_col+8]),
                    float(fields[period_start_col+9]),
                    float(fields[period_start_col+10]),
                    float(fields[period_start_col+11]),
                    float(fields[period_start_col+12]),
                    float(fields[period_start_col+13]),
                    float(fields[period_start_col+14]),
                    float(fields[period_start_col+15]),
                    float(fields[period_start_col+16]),
                    float(fields[period_start_col+17]),
                    float(fields[period_start_col+18]),
                    float(fields[period_start_col+19]),
                    float(fields[period_start_col+20]),
                    float(fields[period_start_col+21]),
                    float(fields[period_start_col+22]),
                    float(fields[period_start_col+23]),
                    float(fields[period_start_col+24]),
                    float(fields[period_start_col+25]),
                    float(fields[period_start_col+26]),
                    float(fields[period_start_col+27]),
                    float(fields[period_start_col+28]),
                    float(fields[period_start_col+29]),
                    float(fields[period_start_col+30]),
                    float(fields[period_start_col+31]),
                    float(fields[period_start_col+32]),
                    float(fields[period_start_col+33]),
                    float(fields[period_start_col+34]),
                    float(fields[period_start_col+35]),
                    float(fields[period_start_col+36]),
                    float(fields[period_start_col+37]),
                    float(fields[period_start_col+38]),
                    float(fields[period_start_col+39]),
                    float(fields[period_start_col+40]),
                    float(fields[period_start_col+41]),
                    float(fields[period_start_col+42]),
                    float(fields[period_start_col+43]),
                    float(fields[period_start_col+44]),
                    float(fields[period_start_col+45]),
                    float(fields[period_start_col+46]),
                    float(fields[period_start_col+47]),
                    float(fields[period_start_col+48]),
                    float(fields[period_start_col+49]),
                    float(fields[period_start_col+50]),
                    float(fields[period_start_col+51]),
                    float(fields[period_start_col+52]),
                    float(fields[period_start_col+53]),
                    float(fields[period_start_col+54]),
                    float(fields[period_start_col+55]),
                    float(fields[period_start_col+56]),
                    float(fields[period_start_col+57]),
                    float(fields[period_start_col+58]),
                    float(fields[period_start_col+59]),
                    float(fields[period_start_col+60]),
                    float(fields[period_start_col+61]),
                    float(fields[period_start_col+62])))
            
            print "Query str:%s" % (cmd)
            cursor.execute(cmd)            
            connection.commit()
    connection.close()
    print "Completed loading event file"
    exit(0)
    
if __name__ == '__main__':
    main()
