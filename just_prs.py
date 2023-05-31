from oakvar import BasePostAggregator
import sqlite3
from sqlite3 import Error
from pathlib import Path
import polars as pl
import urllib

SUM = "sum"
TOTAL = "total"
COUNT = "count"
TITLE = "title"
INVERS = "invers"

class CravatPostAggregator (BasePostAggregator):
    prs:dict = {}
    prs_names:list = []
    # prs5_rsids = []
    sql_get_prs:str = """SELECT name, title, total, invers FROM prs;"""


    def check(self):
        return True


    def setup (self):
        self.sql_file:str = str(Path(__file__).parent) + "/data/prs.sqlite"
        if Path(self.sql_file).exists():
            self.prsconn:sqlite3.Connection = sqlite3.connect(self.sql_file)
            self.prscursor:sqlite3.Cursor = self.prsconn.cursor()
            self.prscursor.execute(self.sql_get_prs)
            rows:tuple = self.prscursor.fetchall()
            for row in rows:
                self.prs_names.append(row[0])
                self.prs[row[0]] = {SUM: 0, COUNT: 0, TITLE: row[1], TOTAL: int(row[2]), INVERS: int(row[3])}

        sql_create:str = """ CREATE TABLE IF NOT EXISTS prs (
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
        self.result_path:Path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.result_conn:sqlite3.Connection = sqlite3.connect(self.result_path)
        self.result_cursor:sqlite3.Cursor = self.result_conn.cursor()
        self.result_cursor.execute(sql_create)
        self.result_conn.commit()
        self.result_cursor.execute("DELETE FROM prs;")


    def get_prs_dataframe(self, name):
        import platform
        sql:str = f"SELECT pos, chrom, effect_allele, weight FROM prs, position, weights WHERE prs.name = '{name}' AND prs.id = weights.prsid AND weights.posid = position.id"
        ol_pl = platform.platform()
        if ol_pl.startswith("Windows"):
            conn_url = f"sqlite://{urllib.parse.quote(self.sql_file)}"
        else:
            conn_url = f"sqlite://{self.sql_file}"
        return pl.read_sql(sql, conn_url)


    def calculate_prs(self, data_df, name):
        prs_df:pl.DataFrame = self.get_prs_dataframe(name)
        prs_df = prs_df.with_column((pl.col('chrom') + pl.col('pos')).alias("key"))
        unite:pl.DataFrame = data_df.join(prs_df, left_on='key', right_on="key")
        unite1 = unite.filter(pl.col("A") == pl.col("effect_allele"))
        unite2 = unite.filter(pl.col("B") == pl.col("effect_allele"))
        res1:pl.Series = unite1.select(pl.col("weight")).sum()
        res2:pl.Series = unite2.select(pl.col("weight")).sum()

        #TODO: negative join for ref homo zygot
        # anti_unite = prs_df.join(data_df, left_on="rsid", right_on="dbsnp__rsid", how="anti")
        # res3 = anti_unite.filter(pl.col("effect_allele") == pl.col("ref")).select(pl.col("weight")).sum() * 2
        return float(res1.item()) + float(res2.item()), unite.shape[0]


    def process_file(self):
        self._close_db_connection()
        data_df = self.get_df("variant", None, 0)
        data_df = data_df.select(['base__pos', 'vcfinfo__zygosity', 'base__ref_base', 'base__alt_base', 'base__chrom'])
        data_df = data_df.with_column(pl.col('vcfinfo__zygosity').fill_null("het"))
        data_df = data_df.with_column((pl.col('base__chrom') + pl.col('base__pos')).alias("key"))

        het_zygot = data_df.filter(pl.col('vcfinfo__zygosity') == 'het')
        het_zygot = het_zygot.with_columns([pl.col('base__ref_base').alias("A"), pl.col('base__alt_base').alias("B")])

        hom_zygot = data_df.filter(pl.col('vcfinfo__zygosity') == 'hom')
        hom_zygot = hom_zygot.with_columns([pl.col('base__alt_base').alias("A"), pl.col('base__alt_base').alias("B")])

        data_df = het_zygot.vstack(hom_zygot)
        for name in self.prs_names:
            sum, count = self.calculate_prs(data_df, name)
            self.prs[name][SUM] = sum
            self.prs[name][COUNT] = count
        self._open_db_connection()


    def cleanup (self):
        if self.result_cursor is not None:
            self.result_cursor.close()
        if self.result_conn is not None:
            self.result_conn.commit()
            self.result_conn.close()
        if self.prscursor is not None:
            self.prscursor.close()
        if self.prsconn is not None:
            self.prsconn.close()

        
    # def annotate (self, input_data):
    #     rsid:str = str(input_data['dbsnp__rsid'])
    #     if rsid == '':
    #         return
    #
    #     if not rsid.startswith("rs"):
    #         rsid = 'rs' + rsid
    #     alt:str = input_data['base__alt_base']
    #     ref:str = input_data['base__ref_base']
    #     chrom:str = input_data['base__chrom']
    #
    #     query:str = f"SELECT prs.name, weights.weight, position.effect_allele FROM position, prs, weights WHERE chrom = '{chrom}'" \
    #             f" AND rsid = '{rsid}' AND weights.posid = position.id AND weights.prsid = prs.id"
    #
    #     self.prscursor.execute(query)
    #     rows:tuple = self.prscursor.fetchall()
    #
    #     if len(rows) == 0:
    #         return
    #
    #     zygot:str = input_data['vcfinfo__zygosity']
    #     for name, weight, allele in rows:
    #         if not (allele == alt or (allele == ref and zygot == 'het')):
    #             continue
    #         weight:float = float(weight)
    #         if allele == alt and zygot == 'hom':
    #             weight = 2 * weight
    #
    #         self.prs[name][SUM] += weight
    #         self.prs[name][COUNT] += 1


    def get_percent(self, name:str, value:float) -> float:
        sql_get_percent:str = f"SELECT 'min', percent, max(value) FROM percentiles, prs WHERE percentiles.prs_id = prs.id AND prs.name = '{name}' AND value <= {value} UNION " \
                        f"SELECT 'max', percent, min(value) FROM percentiles, prs WHERE percentiles.prs_id = prs.id AND prs.name = '{name}' AND value >= {value}"
        self.prscursor.execute(sql_get_percent)
        rows:tuple = self.prscursor.fetchall()
        for row in rows:
            if row[0] == 'min':
                min_percent:float = row[1]
                min_value:float = row[2]
            if row[0] == 'max':
                max_percent:float = row[1]
                max_value:float = row[2]

        if min_value is None:
            return max_percent

        if max_value is None:
            return min_percent

        if abs(min_value - value) > abs(max_value - value):
            return max_percent
        else:
            return min_percent


    def postprocess(self):
        sql:str = """ INSERT INTO prs (name, sum, avg, count, title, total, percent, fraction, invers) VALUES (?,?,?,?,?,?,?,?,?);"""
        for name in self.prs_names:
            avg:float = 0
            if self.prs[name][COUNT] > 0:
                avg = self.prs[name][SUM] / (self.prs[name][COUNT] * 2)
            percent:float = self.get_percent(name, self.prs[name][SUM])
            if type(percent) is not float:
                percent = 0.01
            task:tuple = (name, self.prs[name][SUM], avg, self.prs[name][COUNT], self.prs[name][TITLE], self.prs[name][TOTAL], int(percent * 100),
                    percent, self.prs[name][INVERS])
            self.result_cursor.execute(sql, task)