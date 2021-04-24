from config import config
import sqlite3

config = config
conn = sqlite3.connect(config.DB)

conn.execute('''DROP TABLE IF EXISTS GENOMES;''')
conn.execute('''CREATE TABLE GENOMES
         (GENOME     TEXT    NOT NULL,
          SCINAME    TEXT,
          COMMONNAME TEXT,
          DESCRIPTION TEXT);''')

conn.execute('''DROP TABLE IF EXISTS PROJECTS;''')
conn.execute('''CREATE TABLE PROJECTS
         (PROJECTID    TEXT    NOT NULL,
          OSOURCE    TEXT    NOT NULL,
          DSOURCE    TEXT    NOT NULL,
          SHORTNAME  TEXT    NOT NULL,
          LONGNAME   TEXT    NOT NULL,
          INFO       TEXT    NOT NULL);''')

conn.execute('''DROP TABLE IF EXISTS INDICES;''')
conn.execute('''CREATE TABLE INDICES
         (INDEXID   TEXT    NOT NULL,
          ITER      INT     NOT NULL,
          DATE      NUMERIC NOT NULL,
          PROJECTID TEXT    NOT NULL,
          DSOURCE   TEXT    NOT NULL,
          GENOME    TEXT    NOT NULL,
          FULL      BOOL    NOT NULL,
          SIZE      INT     NOT NULL);''')

conn.execute('''DROP TABLE IF EXISTS FILES;''')
conn.execute('''CREATE TABLE FILES
         (FILEID    TEXT     NOT NULL,
          SIZE      INT      NOT NULL,
          GENOME    TEXT     NOT NULL,
          PROJECTID TEXT     NOT NULL,
          INDEXID   TEXT     NOT NULL,
          SHORTNAME TEXT,
          LONGNAME  TEXT,
          SHORTINFO TEXT,
          LONGINFO  TEXT);''')

print("Database created. Tables initialized PROJECTS, FILES and INDICES")
conn.commit()
conn.close()
