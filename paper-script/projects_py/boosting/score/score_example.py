# how to deal with phe_female

# ng
#phe_config={
#    'phe_female':[],
#}


# ng
#def extract_score(score,phe,cv):
#    if phe in phe_config['phe_female']:
#        pass
#    else:
#        pass

def extract_score(score,phe,cv,sex='female'):
    if sex=='female':
        pass
    else:
        pass


def allstat(phes=[],phes_female=[]):
    for phe in phes:
        if phe in phes_female:
            extract_score(score,phe,sex='female')


# ng
#if __name__=='__main__':
#    phe_config['phe_female']=['bc']