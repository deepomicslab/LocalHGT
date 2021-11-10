from simulation import Parameters

class Batch():

    def __init__(self):
        Parameters.__init__(self)
        self.batch_num = 1
        self.workdir = "/mnt/d/breakpoints/HGT/"
        self.fq_dir = ''
        self.result_dir = ''
        self.localHGT = "/mnt/d/breakpoints/script/pipeline.sh"
        self.hit = 0.1
        self.perfect_hit = 0.03
        self.fq1 = ''
        self.fq2 = ''

    def get_fq_dir(self, fq_dir):
        self.fq_dir = fq_dir

    def get_result_dir(self, result_dir):
        self.result_dir = result_dir

    def get_fq(self):
        self.fq1 = self.fq_dir + '/%s.1.fq'%(sample)
        self.fq2 = self.fq_dir + '/%s.2.fq'%(sample)




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
    for snp_rate in ba.snp_level:
        ba.change_snp_rate(snp_rate)
        for index in range(iteration_times):
            ba.get_ID(index)
            ba.get_fq()
            order = "bash %s %s %s %s %s %s %s %s"%(ba.localHGT, ba.origin_ref, \
            ba.fq1, ba.fq2, ba.sample, ba.result_dir, ba.hit, ba.perfect_hit)
            f = open('%s/work_%s.sh'%(ba.workdir, i%(ba.batch_num)), 'a')
            print (order, file = f)
            f.close()
            i += 1

if __name__ == "__main__":
    batch_snp()