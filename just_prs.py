from oakvar import BasePostAggregator
import sqlite3
from sqlite3 import Error
from pathlib import Path

SUM = "sum"
TOTAL = "total"
COUNT = "count"
TITLE = "title"
INVERS = "invers"

class CravatPostAggregator (BasePostAggregator):
    prs = {}
    prs_names = []
    # prs5_rsids = []
    sql_get_prs = """SELECT name, title, total, invers FROM prs;"""


    def check(self):
        return True


    def setup (self):
        sql_file = str(Path(__file__).parent) + "/data/prs.sqlite"
        if Path(sql_file).exists():
            self.prsconn = sqlite3.connect(sql_file)
            self.prscursor = self.prsconn.cursor()
            self.prscursor.execute(self.sql_get_prs)
            rows = self.prscursor.fetchall()
            for row in rows:
                self.prs_names.append(row[0])
                self.prs[row[0]] = {SUM: 0, COUNT: 0, TITLE: row[1], TOTAL: int(row[2]), INVERS: int(row[3])}

        sql_create = """ CREATE TABLE IF NOT EXISTS prs (
                    id integer NOT NULL PRIMARY KEY,
                    name text,
                    sum float,
                    avg float,
                    count int,
                    title text,
                    total int,
                    percent int,
                    fraction float,
                    invers text
                    )"""
        self.result_path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.longevity_conn = sqlite3.connect(self.result_path)
        self.longevity_cursor = self.longevity_conn.cursor()
        self.longevity_cursor.execute(sql_create)
        self.longevity_conn.commit()
        self.longevity_cursor.execute("DELETE FROM prs;")


    def cleanup (self):
        if self.longevity_cursor is not None:
            self.longevity_cursor.close()
        if self.longevity_conn is not None:
            self.longevity_conn.commit()
            self.longevity_conn.close()
        if self.prscursor is not None:
            self.prscursor.close()
        if self.prsconn is not None:
            self.prsconn.close()

        
    def annotate (self, input_data):
        rsid = str(input_data['dbsnp__rsid'])
        if rsid == '':
            return

        if not rsid.startswith("rs"):
            rsid = 'rs' + rsid
        alt = input_data['base__alt_base']
        ref = input_data['base__ref_base']
        chrom = input_data['base__chrom']

        query = f"SELECT prs.name, weights.weight, position.effect_allele FROM position, prs, weights WHERE chrom = '{chrom}'" \
                f" AND rsid = '{rsid}' AND weights.posid = position.id AND weights.prsid = prs.id"

        self.prscursor.execute(query)
        rows = self.prscursor.fetchall()

        if len(rows) == 0:
            return

        zygot = input_data['vcfinfo__zygosity']
        for name, weight, allele in rows:
            if not (allele == alt or (allele == ref and zygot == 'het')):
                continue
            weight = float(weight)
            if allele == alt and zygot == 'hom':
                weight = 2 * weight

            self.prs[name][SUM] += weight
            self.prs[name][COUNT] += 1
        return {"col1":""}


    def get_percent(self, name, value):
        sql_get_percent = f"SELECT 'min', percent, max(value) FROM percentiles, prs WHERE percentiles.prs_id = prs.id AND prs.name = '{name}' AND value <= {value} UNION " \
                        f"SELECT 'max', percent, min(value) FROM percentiles, prs WHERE percentiles.prs_id = prs.id AND prs.name = '{name}' AND value >= {value}"
        self.prscursor.execute(sql_get_percent)
        rows = self.prscursor.fetchall()
        for row in rows:
            if row[0] == 'min':
                min_percent = row[1]
                min_value = row[2]
            if row[0] == 'max':
                max_percent = row[1]
                max_value = row[2]

        if min_value is None:
            return max_percent

        if max_value is None:
            return min_percent

        if abs(min_value - value) > abs(max_value - value):
            return max_percent
        else:
            return min_percent


    def postprocess(self):
        sql = """ INSERT INTO prs (name, sum, avg, count, title, total, percent, fraction, invers) VALUES (?,?,?,?,?,?,?,?,?);"""
        for name in self.prs_names:
            avg = 0
            if self.prs[name][COUNT] > 0:
                avg = self.prs[name][SUM] / (self.prs[name][COUNT] * 2)
            percent = self.get_percent(name, self.prs[name][SUM])
            if type(percent) is not float:
                percent = 0.01
            task = (name, self.prs[name][SUM], avg, self.prs[name][COUNT], self.prs[name][TITLE], self.prs[name][TOTAL], int(percent * 100),
                    percent, self.prs[name][INVERS])
            self.longevity_cursor.execute(sql, task)