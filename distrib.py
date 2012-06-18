import random 

z_table = {}

def __loadZtable():
    for i in open("z_table","r"):
        i = i.rsplit()
        z_table[float(i[0])] = float(i[1])

def __getZvalue():
    rnd =  random.random()/2
    diffs = [abs(val-rnd) for val in sorted(z_table.values())]
    min_index = diffs.index(min(diffs))
    z_abs = sorted(z_table.keys())[min_index]
    z_value = z_abs * random.choice([-1,1])
    return z_value

def __getValueNorm(var_coef, gene_mean, z_value):
    s = (var_coef * gene_mean)
    xi = (z_value * s) + gene_mean
    return xi

def createNormalData():
    __loadZtable()
    var_coef = 0.2
    gene_mean = 100.001
    for i in xrange(0,300000):
        print __getValueNorm(var_coef, gene_mean, __getZvalue())

def valueNormal(var_coef, gene_mean):
    __loadZtable()
    return __getValueNorm(var_coef, gene_mean, __getZvalue())

