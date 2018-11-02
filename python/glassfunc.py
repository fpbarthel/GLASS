"""
GLASS helper functions
"""

import os, fnmatch
from configparser import ConfigParser

def touch_file(fname, mode=0o666, dir_fd=None, **kwargs):
    """
    Touch function taken from stackoverflow
    Link: https://stackoverflow.com/questions/1158076/implement-touch-using-python
    """
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd, **kwargs)

def build_dict(seq, key):
    """
    Turn an unnamed list of dicts into a nammed list of dicts
    Taken from stackoverflow
    https://stackoverflow.com/questions/4391697/find-the-index-of-a-dict-within-a-list-by-matching-the-dicts-value
    """
    return dict((d[key], dict(d, index=index)) for (index, d) in enumerate(seq))

def dbconfig(filename, section):
    """
    Loads db connection settings in text file
    From http://www.postgresqltutorial.com/postgresql-python/connect/
    """
    # create a parser
    parser = ConfigParser()
    # read config file
    parser.read(filename)
 
    # get section, default to postgresql
    db = {}
    if parser.has_section(section):
        params = parser.items(section)
        for param in params:
            db[param[0]] = param[1]
    else:
        raise Exception('Section {0} not found in the {1} file'.format(section, filename))
 
    return db

def locate(pattern, root = os.curdir):
    """
    Locate all files matching supplied filename pattern in and below
    supplied root directory.
    Taken from: http://code.activestate.com/recipes/499305-locating-files-throughout-a-directory-tree/
    """
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)

## END ##
