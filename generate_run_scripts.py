#!/usr/bin/env python3

from simulation import Parameters

class Batch(Parameters):

    def __init__(self):
        Parameters.__init__(self)
        self.batch_num = 1
        self.workdir = "/mnt/d/breakpoints/HGT/"
        self.fq_dir = ''
        self.result_dir = ''
        self.localHGT = "/mnt/d/breakpoints/script/pipeline.sh"
        self.hit = 0.1
        self.perfect_hit = 0.08
        self.fq1 = ''
        self.fq2 = ''
        self.sample_fasta = ''
        self.LEMON = "/mnt/d/breakpoints/lemon/pipeline.sh"
        self.lemon_outdir = "/mnt/e/HGT/lemon_snp/"

    def get_fq_dir(self, fq_dir):
        self.fq_dir = fq_dir

    def get_result_dir(self, result_dir):
        self.result_dir = result_dir

    def get_fq(self):
        self.fq1 = self.fq_dir + '/%s.1.fq'%(self.sample)
        self.fq2 = self.fq_dir + '/%s.2.fq'%(self.sample)

    def get_minor_order(self):
        order = "/usr/bin/time -v -o %s/%s.time bash %s %s/%s.fa %s %s %s %s %s %s"%(self.result_dir, \
        self.sample, self.localHGT, self.fq_dir, self.sample, \
        self.fq1, self.fq2, self.sample, self.result_dir, self.hit, self.perfect_hit)
        return order

    def get_normal_order(self):
        order = "/usr/bin/time -v -o %s/%s.time bash %s %s %s %s %s %s %s %s"%(self.result_dir, self.sample,\
         self.localHGT, self.origin_ref, \
        self.fq1, self.fq2, self.sample, self.result_dir, self.hit, self.perfect_hit)
        return order 

    def get_lemon_order(self):
        order = f"/usr/bin/time -v -o {self.lemon_outdir}/{self.sample}.time bash {self.LEMON} {self.origin_ref}\
         {self.fq1} {self.fq2} {self.sample} {self.lemon_outdir}"
        return order

    def change_ID(self, new_id):
        self.sample = new_id


def batch_snp():
    ba = Batch()
    ba.get_fq_dir("/mnt/d/breakpoints/HGT/uhgg_snp/")
    ba.get_result_dir("/mnt/d/breakpoints/HGT/uhgg_snp_results/")

    o = open(ba.workdir + 'work.sh', 'w')
    for h in range(ba.batch_num):
        f = open(ba.workdir +'work_%s.sh'%(h), 'w')
        f.close()
        print ('nohup sh %s/work_%s.sh &>%s/log%s&'%(ba.workdir, h, ba.workdir, h), file = o)

    i = 1
    index = 0
    for snp_rate in ba.snp_level[1:-1]:
    # for snp_rate in [0.08, 0.09]:
        ba.change_snp_rate(snp_rate)
        for index in range(ba.iteration_times):
            ba.get_ID(index)
            ba.get_fq()
            order = ba.get_minor_order()
            f = open('%s/work_%s.sh'%(ba.workdir, i%(ba.batch_num)), 'a')
            print (order, file = f)
            f.close()
            i += 1

def batch_cami():
    ba = Batch()
    ba.get_fq_dir("/mnt/d/breakpoints/HGT/uhgg_snp/")
    ba.get_result_dir("/mnt/d/breakpoints/HGT/uhgg_snp_results/")

    f = open("/mnt/d/breakpoints/HGT/run_lemon_cami.sh", 'w')
    h = open("/mnt/d/breakpoints/HGT/run_localHGT_cami.sh", 'w')

    i = 1
    index = 0
    # for snp_rate in [0.09, 0.07, 0.05, 0.03, 0.01]:
    for snp_rate in [0.01, 0.02, 0.03, 0.04, 0.05]:
        ba.change_snp_rate(snp_rate)
        index = 0
        ba.get_ID(index)
        for level in ba.complexity_level:
            cami_ID = ba.sample + '_' + level
            ba.change_ID(cami_ID)
            ba.get_fq()
            order = ba.get_lemon_order()
            print (order, file = f)
            order = ba.get_normal_order()
            print (order, file = h)
            ba.get_ID(index) # refresh ID

    f.close()
    h.close()


if __name__ == "__main__":
    # batch_snp()
    batch_cami()