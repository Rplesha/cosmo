from __future__ import print_function, absolute_import, division

import os

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import ForeignKey, Column, Index, Integer, String, Float, Boolean, Numeric, BigInteger, Text
from sqlalchemy.dialects import mysql
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker, relationship, backref

try:
    import yaml
except ImportError:
    from .yaml import yaml

__all__ = ['open_settings', 'load_connection']

Base = declarative_base()

#-------------------------------------------------------------------------------

def open_settings(config_file=None):
    """ Parse config file and load settings

    If no config file is supplied, the configuration file will assume to be
    located at '~/configure.yaml'.

    Parameters
    ----------
    config_file : str, optional
        yaml file containing configuration settings.

    Returns
    -------
    settings : dict
        dictionary of all settings

    """

    config_file = config_file or os.path.join(os.environ['HOME'], "configure.yaml")

    with open(config_file, 'r') as f:
        settings = yaml.load(f)

    return settings

#-------------------------------------------------------------------------------

def load_connection(connection_string, echo=False):
    """Create and return a connection to the database given in the
    connection string.

    Parameters
    ----------
    connection_string : str
        A string that points to the database conenction.  The
        connection string is in the following form:
        dialect+driver://username:password@host:port/database
    echo : bool
        Show all SQL produced.

    Returns
    -------
    session : sesson object
        Provides a holding zone for all objects loaded or associated
        with the database.
    engine : engine object
        Provides a source of database connectivity and behavior.
    """

    engine = create_engine(connection_string, echo=echo)
    Session = sessionmaker(bind=engine)

    return Session, engine

#-------------------------------------------------------------------------------

class Darks(Base):
    __tablename__ = "darks"

    id = Column(Integer, primary_key=True)

    obsname = Column(String(30))
    rootname = Column(String(9))
    detector = Column(String(4))
    date = Column(Float)
    dark = Column(Float)
    ta_dark = Column(Float)
    latitude = Column(Float)
    longitude = Column(Float)
    sun_lat = Column(Float)
    sun_lon = Column(Float)
    temp = Column(Float)

    file_id = Column(Integer, ForeignKey('files.id'))
    #file = relationship("Files", backref=backref('lampflash', order_by=id))

#-------------------------------------------------------------------------------

class Files(Base):
    __tablename__ = 'files'

    id = Column(Integer, primary_key=True)

    path = Column(String(70))
    name = Column(String(40))
    rootname = Column(String(9))

    __table_args__ = (Index('idx_fullpath', 'path', 'name', unique=True), )
    __table_args__ = (Index('idx_rootname', 'rootname'), )

#-------------------------------------------------------------------------------

class Lampflash(Base):
    __tablename__ = 'lampflash'

    id = Column(Integer, primary_key=True)

    date = Column(Float)
    rootname = Column(String(9))
    proposid = Column(Integer)
    detector = Column(String(4))
    segment = Column(String(4))
    opt_elem = Column(String(7))
    cenwave = Column(Integer)
    fppos = Column(Integer)
    lamptab = Column(String(30))
    flash = Column(Integer)
    x_shift = Column(Float)
    y_shift = Column(Float)
    found = Column(Boolean)

    file_id = Column(Integer, ForeignKey('files.id'))
    __table_args__ = (Index('idx_rootname', 'rootname', unique=False), )
    #file = relationship("Files", backref=backref('lampflash', order_by=id))

#-------------------------------------------------------------------------------

class Headers(Base):
    __tablename__ = "headers"

    id = Column(Integer, primary_key=True)
    filetype = Column(String(67))
    instrume = Column(String(3))
    rootname = Column(String(9))
    imagetyp = Column(String(20))
    targname = Column(String(67))
    ra_targ = Column(Float(20))
    dec_targ = Column(Float(20))
    proposid = Column(Integer)
    qualcom1 = Column(String(67))
    qualcom2 = Column(String(67))
    qualcom3 = Column(String(67))
    quality = Column(String(67))
    opus_ver = Column(String(30))
    postarg1 = Column(Float(32))
    postarg2 = Column(Float)
    cal_ver = Column(String(30))
    proctime = Column(Float)

    obstype = Column(String(20))
    obsmode = Column(String(20))
    exptype = Column(String(20))
    detector = Column(String(20))
    segment = Column(String(20))
    detecthv = Column(String(20))
    life_adj = Column(Integer)
    fppos = Column(Integer)
    exp_num = Column(Integer)
    cenwave = Column(Integer)
    propaper = Column(String(20))
    apmpos = Column(String(20))
    aperxpos = Column(Float)
    aperypos = Column(Float)
    aperture = Column(String(4))
    opt_elem = Column(String(7))
    shutter = Column(String(20))
    extended = Column(String(20))
    obset_id = Column(String(2))
    asn_id = Column(String(9))
    asn_tab = Column(String(18))
    randseed = Column(BigInteger)
    asn_mtyp = Column(String(20))
    overflow = Column(Integer)
    nevents = Column(Integer)
    neventsa = Column(Float)
    neventsb = Column(Float)
    dethvla = Column(Integer)
    dethvlb = Column(Integer)
    deventa = Column(Float)
    deventb = Column(Float)
    feventa = Column(Float)
    feventb = Column(Float)
    hvlevela = Column(Integer)
    hvlevelb = Column(Integer)
    dpixel1a = Column(Float)
    dpixel1b = Column(Float)
    date_obs = Column(String(10))
    time_obs = Column(String(8))
    #expstart = Column(Numeric(8, 3))
    #expend = Column(Numeric(8, 3))
    expstart = Column(Float)
    expend = Column(Float)
    exptime = Column(Float)
    numflash = Column(Integer)
    ra_aper = Column(Float)
    dec_aper = Column(Float)
    shift1a = Column(Float)
    shift2a = Column(Float)
    shift1b = Column(Float)
    shift2b = Column(Float)
    shift1c = Column(Float)
    shift2c = Column(Float)

    sp_loc_a = Column(Float)
    sp_loc_b = Column(Float)
    sp_loc_c = Column(Float)
    sp_nom_a = Column(Float)
    sp_nom_b = Column(Float)
    sp_nom_c = Column(Float)
    sp_off_a = Column(Float)
    sp_off_b = Column(Float)
    sp_off_c = Column(Float)
    sp_err_a = Column(Float)
    sp_err_b = Column(Float)
    sp_err_c = Column(Float)

    #NUV keywords
    dethvl = Column(Float)

    file_id = Column(Integer, ForeignKey('files.id'))
    #file = relationship("Files", backref=backref('headers', order_by=id))

    __table_args__ = (Index('idx_rootname', 'rootname', unique=False), )
    __table_args__ = (Index('idx_config', 'segment', 'fppos', 'cenwave', 'opt_elem', unique=False), )

#-------------------------------------------------------------------------------

class Data(Base):
    __tablename__ = "data"

    id = Column(Integer, primary_key=True)

    flux_mean = Column(Float)
    flux_max = Column(Float)
    flux_std = Column(Float)
    wl_min = Column(Float)
    wl_max = Column(Float)


    file_id = Column(Integer, ForeignKey('files.id'))
    #file = relationship("Files", backref=backref('Data', order_by=id))

#-------------------------------------------------------------------------------

class Stims(Base):
    """Record location of all STIM pulses"""
    __tablename__ = "stims"

    id = Column(Integer, primary_key=True)

    time = Column(Float)
    rootname = Column(String(9))
    abs_time = Column(Float)
    stim1_x = Column(Float)
    stim1_y = Column(Float)
    stim2_x = Column(Float)
    stim2_y = Column(Float)
    counts = Column(Float)
    segment = Column(String(4))
    file_id = Column(Integer, ForeignKey('files.id'))

    __table_args__ = (Index('idx_rootname', 'rootname', unique=False), )
    #file = relationship("Files", backref=backref('Stims', order_by=id))

#-------------------------------------------------------------------------------

class Phd(Base):
    __tablename__ = 'phd'

    id = Column(Integer, primary_key=True)

    pha_0 = Column(Integer)
    pha_1 = Column(Integer)
    pha_2 = Column(Integer)
    pha_3 = Column(Integer)
    pha_4 = Column(Integer)
    pha_5 = Column(Integer)
    pha_6 = Column(Integer)
    pha_7 = Column(Integer)
    pha_8 = Column(Integer)
    pha_9 = Column(Integer)
    pha_10 = Column(Integer)
    pha_11 = Column(Integer)
    pha_12 = Column(Integer)
    pha_13 = Column(Integer)
    pha_14 = Column(Integer)
    pha_15 = Column(Integer)
    pha_16 = Column(Integer)
    pha_17 = Column(Integer)
    pha_18 = Column(Integer)
    pha_19 = Column(Integer)
    pha_20 = Column(Integer)
    pha_21 = Column(Integer)
    pha_22 = Column(Integer)
    pha_23 = Column(Integer)
    pha_24 = Column(Integer)
    pha_26 = Column(Integer)
    pha_25 = Column(Integer)
    pha_27 = Column(Integer)
    pha_28 = Column(Integer)
    pha_29 = Column(Integer)
    pha_30 = Column(Integer)
    pha_31 = Column(Integer)

    file_id = Column(Integer, ForeignKey('files.id'))
    #file = relationship("Files", backref=backref('Phd', order_by=id))

#-------------------------------------------------------------------------------

class Gain(Base):
    __tablename__ = 'gain'

    id = Column(BigInteger, primary_key=True)

    x = Column(Integer)
    y = Column(Integer)
    gain = Column(Float)
    counts = Column(Float)
    std = Column(Float)
    segment = Column(String(4))
    dethv = Column(Integer)
    expstart = Column(Float)

    file_id = Column(Integer, ForeignKey('files.id'))
    __table_args__ = (Index('coord', 'x', 'y', unique=False), )
    #file = relationship("Files", backref=backref('Gain', order_by=id))

#-------------------------------------------------------------------------------

class sptkeys(Base):
    __tablename__ = 'spt'

    id = Column(Integer, primary_key=True)

    #spt file keywords
    rootname = Column(String(9))
    proc_typ = Column(String(20)) # primary extention
    lomfstp = Column(Float) #2 ext, focus in spreadsheet
    lapxlvdt = Column(Integer) #2 ext, aper_disp in spreadsheet
    lapdlvdt = Column(Integer) #2 ext, aper_xdisp in spreadsheet
    lom1posc = Column(Integer) #2 ext, osm1_coarse in spreadsheet
    lom2posc = Column(Integer) #2 ext, osm2_coarse in spreadsheet
    lom1posf = Column(Integer) #2 ext, osm1_fine in spreadsheet
    lom2posf = Column(Integer) #2 ext, osm2_fine in spreadsheet
    ldcampat = Column(Float)
    ldcampbt = Column(Float)
    lmmcetmp = Column(Float)

    file_id = Column(Integer, ForeignKey('files.id'))

    __table_args__ = (Index('idx_rootname', 'rootname', unique=False), )

#-------------------------------------------------------------------------------
