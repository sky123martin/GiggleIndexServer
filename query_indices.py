import subprocess
from config import config
import sqlite3
from sys import argv

config = config
conn = sqlite3.connect(config.DB)

# def query_interval(intv):
#     giggle search -i data/UCSC_rn6.1.d -r 1:200457776-200457776

# def query_interval(intv):