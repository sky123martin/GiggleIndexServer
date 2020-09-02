from config import config
import sqlite3

config = config
conn = sqlite3.connect(config.DB)

def update_genome(source, genome):
    pass


if __name__ == "__main__":
    out = conn.execute("SELECT SOURCE, GENOME from FILES")
    for row in out:
        source = row[0]
        genome = row[1]
        update_genome(source, genome)