#!/usr/bin/env python
"""
Name: create_bbp_db.py
Purpose: create the db-schema for bbp
The schema will be defined here, combined with code to create the tables.
Later enhancmements will be to export from sqlite and recreate from export format

Table naming conventions are:
- SQL CMDS all capitals
- Table Names (Noun in plural form)
- Row Names (First letter capital with underscores

Author: Philip Maechling
Created: Jul 21, 2012
$Id: make_and_load_bbp_xml.py 1625 2016-04-08 18:41:45Z fsilva $
"""

# Import Python modules
import os
import sys
import time
import sqlite3
from optparse import OptionParser
from xml.dom.minidom import parseString

def create_bbp_db(sqlitedb):
    """
    This function creates the BBP DB database
    """
    if os.path.exists(sqlitedb):
        print 'Warning...This script will delete existing BBP DB...'
        print "Press <Enter> to continue"
        try:
            os.read(1, 1)
        except KeyboardInterrupt:
            print "Cancelled!"
            sys.exit(0)

    print "Starting BBP_DB table creation..."

    connection = sqlite3.connect(sqlitedb)
    cursor = connection.cursor()

    # Drop all the tables so they do not exist when we try to create them
    strcmd = '\
    DROP TABLE IF EXISTS Versions'
    cursor.execute(strcmd)
    strcmd = '\
    DROP TABLE IF EXISTS Studys'
    cursor.execute(strcmd)
    strcmd = '\
    DROP TABLE IF EXISTS Modules'
    cursor.execute(strcmd)
    strcmd = '\
    DROP TABLE IF EXISTS Velocity_Models'
    cursor.execute(strcmd)
    strcmd = '\
    DROP TABLE IF EXISTS Events'
    cursor.execute(strcmd)
    strcmd = '\
    DROP TABLE IF EXISTS RSNs'
    cursor.execute(strcmd)
    connection.commit()

    # -------- Create Versions Table --------
    strcmd = '\
    CREATE TABLE IF NOT EXISTS Versions \
    (Schema_ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, \
    Schema_Version TEXT NOT NULL, \
    Schema_Date NUMERIC)'

    cursor.execute(strcmd)
    connection.commit()

    # -------- Create Studys Table --------
    # define a short name as primary key
    strcmd = '\
    CREATE TABLE IF NOT EXISTS Studys \
    (Study_Name TEXT NOT NULL PRIMARY KEY, \
    Study_Description TEXT NOT NULL, \
    Study_Contact TEXT NOT NULL)'

    cursor.execute(strcmd)
    connection.commit()

    # -------- Create Modules Table ---------
    strcmd = '\
    CREATE TABLE IF NOT EXISTS Modules \
    (Module_ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, \
    Module_Name TEXT NOT NULL UNIQUE, \
    Module_Descript TEXT)'

    cursor.execute(strcmd)
    connection.commit()

    # -------- Create Velocity Models Table ---------
    strcmd = '\
    CREATE TABLE IF NOT EXISTS Velocity_Models \
    (Vel_Model_ID INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT, \
    Vel_Model_Name TEXT NOT NULL UNIQUE, \
    Vel_Model_Description TEXT)'

    cursor.execute(strcmd)
    connection.commit()

    # -------- Create Event Table ---------
    strcmd = '\
    CREATE TABLE IF NOT EXISTS Events \
    (EQID INTEGER NOT NULL PRIMARY KEY, \
    Event_Name TEXT NOT NULL, \
    eq_year INTEGER, \
    eq_mody INTEGER, \
    eq_hrmn INTEGER, \
    eq_mag REAL, \
    hypo_lat REAL , \
    hypo_lon REAL , \
    hypo_depth REAL )'

    cursor.execute(strcmd)
    connection.commit()

    # ---- Create RSN Table with each Station linked to a specific event ----
    strcmd = '\
    CREATE TABLE IF NOT EXISTS RSNs \
    (RECID INTEGER NOT NULL, \
    EQID INTEGER NOT NULL, \
    StationID TEXT NOT NULL, \
    Study_Name  TEXT NOT NULL, \
    SLAT  REAL, \
    SLONG REAL, \
    Rrup REAL, \
    VS30 REAL, \
    Lowpass_Freq REAL, \
    Highpass_Freq REAL, \
    n_file TEXT NOT NULL UNIQUE, \
    e_file TEXT NOT NULL UNIQUE, \
    z_file TEXT NOT NULL UNIQUE, \
    PGAr REAL, \
    p0 REAL, \
    p1 REAL, \
    p2 REAL, \
    p3 REAL, \
    p4 REAL, \
    p5 REAL, \
    p6 REAL, \
    p7 REAL, \
    p8 REAL, \
    p9 REAL, \
    p10 REAL, \
    p11 REAL, \
    p12 REAL, \
    p13 REAL, \
    p14 REAL, \
    p15 REAL, \
    p16 REAL, \
    p17 REAL, \
    p18 REAL, \
    p19 REAL, \
    p20 REAL, \
    p21 REAL, \
    p22 REAL, \
    p23 REAL, \
    p24 REAL, \
    p25 REAL, \
    p26 REAL, \
    p27 REAL, \
    p28 REAL, \
    p29 REAL, \
    p30 REAL, \
    p31 REAL, \
    p32 REAL, \
    p33 REAL, \
    p34 REAL, \
    p35 REAL, \
    p36 REAL, \
    p37 REAL, \
    p38 REAL, \
    p39 REAL, \
    p40 REAL, \
    p41 REAL, \
    p42 REAL, \
    p43 REAL, \
    p44 REAL, \
    p45 REAL, \
    p46 REAL, \
    p47 REAL, \
    p48 REAL, \
    p49 REAL, \
    p50 REAL, \
    p51 REAL, \
    p52 REAL, \
    p53 REAL, \
    p54 REAL, \
    p55 REAL, \
    p56 REAL, \
    p57 REAL, \
    p58 REAL, \
    p59 REAL, \
    p60 REAL, \
    p61 REAL, \
    p62 REAL, \
    PRIMARY KEY (RECID,EQID), \
    FOREIGN KEY (EQID) REFERENCES Events(EQID))'
    #print str
    cursor.execute(strcmd)
    connection.commit()

    #-- Wrapup --
    cursor.close()
    connection.close()
    print "Created BBP_DB Tables"

def load_bbp_xml(sqlitedb, dbfile):
    """
    This function loads the schema from dbfile into the sqlitedb
    """
    # open the xml file for reading:
    file = open(dbfile, 'r')
    # convert to string:
    data = file.read()
    # close file because we dont need it anymore:
    file.close()
    # parse the xml you got from the file
    dom = parseString(data)

    #
    # Connect to sqlite3db
    #
    connection = sqlite3.connect(sqlitedb)
    cursor = connection.cursor()

    #
    # Load "Versions" Table
    #
    # retrieve the first xml tag (<tag>data</tag>) that the parser
    # finds with name tagName:
    xmltag = dom.getElementsByTagName('version')[0].toxml()
    # strip off the tag (<tag>data</tag>  --->   data):
    xmldata = xmltag.replace('<version>', '').replace('</version>', '')
    vstring = xmldata

    # Populate the Versions Table
    today = time.strftime("%A, %B %d, %Y")
    cursor.execute('INSERT INTO Versions VALUES(null, ?, ?)', (vstring, today))
    connection.commit()
    print "Loaded Version Table value: %s" % (vstring)

    #
    # Load "Studys" Table
    #
    # retrieve the first xml tag (<tag>data</tag>) that the parser
    # finds with name tagName:
    studylist = dom.getElementsByTagName('study')
    for s in studylist:
        snametag = s.getElementsByTagName('study_name')[0].toxml()
        sname_data = snametag.replace('<study_name>', '').replace('</study_name>', '')
        sdescripttag = s.getElementsByTagName('study_descript')[0].toxml()
        sdescript_data = sdescripttag.replace('<study_descript>', '').replace('</study_descript>', '')
        study_contacttag = s.getElementsByTagName('study_contact')[0].toxml()
        study_contact_data = study_contacttag.replace('<study_contact>', '').replace('</study_contact>', '')

        # Populate the Studys Table
        cursor.execute('INSERT INTO Studys VALUES(?, ?, ?)', (sname_data,
                                                              sdescript_data,
                                                              study_contact_data))
        connection.commit()
        print "Loaded Studys Table value: %s" % (sname_data)

    #
    # Load "Modules" Table
    #
    # retrieve the first xml tag (<tag>data</tag>) that the parser
    # finds with name tagName:
    studylist = dom.getElementsByTagName('module')
    for s in studylist:
        snametag = s.getElementsByTagName('module_name')[0].toxml()
        sname_data = snametag.replace('<module_name>', '').replace('</module_name>', '')
        sdescripttag = s.getElementsByTagName('module_descript')[0].toxml()
        sdescript_data = sdescripttag.replace('<module_descript>', '').replace('</module_descript>', '')
        cursor.execute('INSERT INTO Modules VALUES(null, ?, ?)',
                       (sname_data, sdescript_data))
        print "Loaded Modules Table value: %s" % (sname_data)
        connection.commit()

    #
    # Load "Velocity Model" Table
    #
    # retrieve the first xml tag (<tag>data</tag>) that the parser finds
    # with name tagName:
    studylist = dom.getElementsByTagName('velocity_model')
    for s in studylist:
        snametag = s.getElementsByTagName('vmod_name')[0].toxml()
        sname_data = snametag.replace('<vmod_name>', '').replace('</vmod_name>', '')
        sdescripttag = s.getElementsByTagName('vmod_descript')[0].toxml()
        sdescript_data = sdescripttag.replace('<vmod_descript>', '').replace('</vmod_descript>', '')
        cursor.execute('INSERT INTO Velocity_Models VALUES(null, ?, ?)',
                       (sname_data, sdescript_data))
        print "Loaded Velocity Models Table value: %s" % (sname_data)
        connection.commit()

    #
    # Load "Event" Table
    #
    studylist = dom.getElementsByTagName('event')
    for s in studylist:
        snametag = s.getElementsByTagName('eqid')[0].toxml()
        eqid = snametag.replace('<eqid>', '').replace('</eqid>', '')
        sdescripttag = s.getElementsByTagName('eventname')[0].toxml()
        eqname = sdescripttag.replace('<eventname>', '').replace('</eventname>', '')
        cursor.execute('INSERT INTO Events (EQID, Event_Name) VALUES(?,?)',
                       (eqid, eqname))
        connection.commit()
        print "Loaded Event Table value: %s" % (eqname)
        connection.commit()

    # -- Wrapup --
    cursor.close()
    connection.close()
    print "Created BBP_DB tables!"
    exit(0)

def main():
    """
    Parse command-line options and call create and load functions
    """
    sqlitedb = None
    dbfile = None

    usage = "usage: %s [options]" % (sys.argv[0])
    parser = OptionParser(usage)

    parser.add_option("-d", "--db_file", dest="sqlitedb",
                      help="File containing SQLite3 DB")
    parser.add_option("-f", "--db_xml_file", dest="dbfile",
                      help="File containing basic configuration for bbp db")

    #
    # Start of parsing input parameters to db query
    #
    (options, args) = parser.parse_args()

    #
    # Point to user supplied DB is needed
    #
    if options.sqlitedb is not None:
        sqlitedb = options.sqlitedb
        print "Creating BBP DB file: %s" % (sqlitedb)

    if options.dbfile is not None:
        dbfile = options.dbfile
        print "Loading BBP XML file: %s\n" % (dbfile)

    if sqlitedb is None or dbfile is None:
        print "Please provide both database and event files!"
        sys.exit(1)

    create_bbp_db(sqlitedb)
    load_bbp_xml(sqlitedb, dbfile)
    sys.exit(0)

if __name__ == '__main__':
    main()
