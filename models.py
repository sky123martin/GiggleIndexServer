from config import config
import sqlite3

config = config
conn = sqlite3.connect(config.DB)

conn.execute('''DROP TABLE IF EXISTS INDICES;''')
conn.execute('''CREATE TABLE INDICES
         (NAME      TEXT    NOT NULL,
          ITER      INT     NOT NULL,
          DATE      NUMERIC NOT NULL,
          SOURCE    TEXT    NOT NULL,
          GENOME    EXT     NOT NULL,
          FULL      NUMERIC NOT NULL,
          SIZE      INT     NOT NULL);''')

conn.execute('''DROP TABLE IF EXISTS FILES;''')
conn.execute('''CREATE TABLE FILES
         (NAME      TEXT     NOT NULL,
          DATE      NUMERIC  NOT NULL,
          SOURCE    TEXT     NOT NULL,
          GENOME    TEXT     NOT NULL,
          SIZE      INT      NOT NULL,
          INDEXNAME TEXT     NOT NULL,
          SHORTNAME TEXT,
          LONGNAME  TEXT,
          SHORTINFO TEXT,
          LONGINFO  TEXT);''')

print("Database created. Tables initialized FILES and INDICES")
conn.commit()
conn.close()
