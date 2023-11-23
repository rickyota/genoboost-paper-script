import seaborn as sns
from logging import getLogger


logger = getLogger(__name__)


def _methodd_capital():
    methodd = {
        'boosting_integrate': 'GenoBoost',
        'boosting_add': 'Additive GenoBoost',
        'boosting_nonadd': 'Non-additive GenoBoost',
        'boosting_nonadd-effeps2': 'Non-additive GenoBoost adjustment 2',
        'boosting_nonadd-noeffeps': 'Non-additive GenoBoost without limiting $s_2$',
        'boosting_free': 'Free GenoBoost',
        'boosting_integrate-ex-chr6': 'GenoBoost without Chr6',
        'boosting_add-ex-chr6': 'Additive GenoBoost without Chr6',
        'boosting_nonadd-ex-chr6': 'Non-additive GenoBoost without Chr6',
        'boosting_integrate-noneur': 'GenoBoost nonEUR',
        'boosting_add-noneur': 'Additive GenoBoost nonEUR',
        'boosting_nonadd-noneur': 'Non-additive GenoBoost nonEUR',
        'snpboost': 'snpboost',
        'snpboost_noneur': 'snpboost nonEUR',
        'snpboost_ex-chr6': 'snpboost without Chr6',
        'snpnet': 'snpnet',
        'snpnet_nonadd': 'Non-additive snpnet',
        'snpnet_noneur': 'snpnet nonEUR',
        'snpnet_ex-chr6': 'snpnet without Chr6',
        'lassosum': 'lassosum',
        'lassosum_noneur': 'lassosum nonEUR',
        'lassosum_ex-chr6': 'lassosum without Chr6',
        'ldpred': 'LDpred',
        'ldpred_noneur': 'LDpred nonEUR',
        'ldpred_ex-chr6': 'LDpred without Chr6',
        'prscs': 'PRS-CS',
        'prscs_noneur': 'PRS-CS nonEUR',
        'prscs_ex-chr6': 'PRS-CS without Chr6',
        'sbayesr': 'SBayesR',
        'sbayesr_noneur': 'SBayesR nonEUR',
        'sbayesr_ex-chr6': 'SBayesR without Chr6',
        'clump': 'C+T',
        'clump_noneur': 'C+T nonEUR',
        'clump_ex-chr6': 'C+T without Chr6',
    }
    return methodd


def _methodd_capital2():
    methodd = {
        'boosting_integrate': 'GenoBoost',
        'boosting_add': 'Additive\nGenoBoost',
        'boosting_nonadd': 'Non-additive\nGenoBoost',
        'boosting_nonadd-effeps2': 'Non-additive\nGenoBoost\nadjustment 2',
        'boosting_nonadd-noeffeps': 'Non-additive\nGenoBoost\nw/o limiting $s_2$',
        'boosting_free': 'Free\nGenoBoost',
        'boosting_integrate-ex-chr6': 'GenoBoost\nw/o Chr6',
        'boosting_add-ex-chr6': 'Additive\nGenoBoost\nw/o Chr6',
        'boosting_nonadd-ex-chr6': 'Non-additive\nGenoBoost\nw/o Chr6',
        'boosting_integrate-noneur': 'GenoBoost\nnonEUR',
        'boosting_add-noneur': 'Additive\nGenoBoost\nnonEUR',
        'boosting_nonadd-noneur': 'Non-additive\nGenoBoost\nnonEUR',
        'snpboost': 'snpboost',
        'snpboost_noneur': 'snpboost\nnonEUR',
        'snpboost_ex-chr6': 'snpboost\nw/o Chr6',
        'snpnet': 'snpnet\n(LASSO)',
        'snpnet_nonadd': 'Non-additive\nsnpnet',
        'snpnet_noneur': 'snpnet\nnonEUR',
        'snpnet_ex-chr6': 'snpnet\nw/o Chr6',
        'lassosum': 'lassosum',
        'lassosum_noneur': 'lassosum\nnonEUR',
        'lassosum_ex-chr6': 'lassosum\nw/o Chr6',
        'ldpred': 'LDpred',
        'ldpred_noneur': 'LDpred\nnonEUR',
        'ldpred_ex-chr6': 'LDpred\nw/o Chr6',
        'prscs': 'PRS-CS',
        'prscs_noneur': 'PRS-CS\nnonEUR',
        'prscs_ex-chr6': 'PRS-CS\nw/o Chr6',
        'sbayesr': 'SBayesR',
        'sbayesr_noneur': 'SBayesR\nnonEUR',
        'sbayesr_ex-chr6': 'SBayesR\nw/o Chr6',
        'clump': 'C+T',
        'clump_noneur': 'C+T\nnonEUR',
        'clump_ex-chr6': 'C+T\nw/o Chr6',
    }
    return methodd


def _romand_method():
    romand = {
        'boosting_integrate': 'i',
        'boosting_nonadd': 'A',
        'boosting_add': 'B',
        'boosting_nonadd-effeps2': "A2",
        'boosting_nonadd-noeffeps': "A3",
        'boosting_free': 'C',
        'boosting_integrate-ex-chr6': "i'",
        'boosting_nonadd-ex-chr6': "A'",
        'boosting_add-ex-chr6': "B'",
        'boosting_integrate-noneur': "i'",
        'boosting_nonadd-noneur': "A'",
        'boosting_add-noneur': "B'",
        'snpboost': 'ii',
        'snpboost_noneur': "iii'",
        'snpboost_ex-chr6': "iii'",
        'snpnet': 'iii',
        'snpnet_nonadd': "iii'",
        'snpnet_noneur': "iii'",
        'snpnet_ex-chr6': "iii'",
        'lassosum': 'iv',
        'lassosum_noneur': "iv'",
        'lassosum_ex-chr6': "iv'",
        'ldpred': 'v',
        'ldpred_noneur': "v'",
        'ldpred_ex-chr6': "v'",
        'prscs': 'vi',
        'prscs_noneur': "vi'",
        'prscs_ex-chr6': "vi'",
        'sbayesr': 'vii',
        'sbayesr_noneur': "vii'",
        'sbayesr_ex-chr6': "vii'",
        'clump': 'viii',
        'clump_noneur': "viii'",
        'clump_ex-chr6': "viii'",
    }
    return romand


# def _romand_method():
#    romand = {
#        'boosting_integrate': 'i',
#        'boosting_nonadd': 'a',
#        'boosting_add': 'b',
#        'boosting_nonadd-effeps2': "a2",
#        'boosting_nonadd-noeffeps': "a3",
#        'boosting_integrate-ex-chr6': "i'",
#        'boosting_nonadd-ex-chr6': "a'",
#        'boosting_add-ex-chr6': "b'",
#        'boosting_integrate-noneur': "i'",
#        'boosting_nonadd-noneur': "a'",
#        'boosting_add-noneur': "b'",
#        'snpboost': 'xx',
#        'snpnet': 'ii',
#        'snpnet_nonadd': "ii'",
#        'snpnet_noneur': "ii'",
#        'snpnet_ex-chr6': "ii'",
#        'lassosum': 'iii',
#        'lassosum_noneur': "iii'",
#        'lassosum_ex-chr6': "iii'",
#        'ldpred': 'iv',
#        'ldpred_noneur': "iv'",
#        'ldpred_ex-chr6': "iv'",
#        'prscs': 'v',
#        'prscs_noneur': "v'",
#        'prscs_ex-chr6': "v'",
#        'sbayesr': 'vi',
#        'sbayesr_noneur': "vi'",
#        'sbayesr_ex-chr6': "vi'",
#        'clump': 'vii',
#        'clump_noneur': "vii'",
#        'clump_ex-chr6': "vii'",
#    }
#    return romand


def palette_hls_2_methods():
    # create by my self
    # 2 red: boosting
    # 1 orange : lasso (raw genot)
    # 2 green : lassosum, clump
    # 3 blue : bayes (ldpred, sbayesr, prscs)

    cb = sns.color_palette("colorblind", 10)
    # ref: https://qiita.com/SaitoTsutomu/items/c79c9973a92e1e2c77a7
    # ref: sns_palette.ipynb
    hls = sns.hls_palette(24, l=0.5, s=0.8)

    paletted = {
        'boosting_integrate': cb[3],
        'boosting_add': cb[6],
        'boosting_nonadd': cb[1],
        #'boosting_add': cb[1],
        #'boosting_nonadd': cb[6],
        'boosting_nonadd-effeps2': cb[8],
        'boosting_nonadd-noeffeps': cb[7],
        'boosting_free':hls[2],
        # 'boosting_nonadd-noeffeps': cb[6],
        'boosting_integrate-ex-chr6': cb[3],
        'boosting_add-ex-chr6': cb[6],
        'boosting_nonadd-ex-chr6': cb[1],
        'boosting_integrate-noneur': cb[3],
        'boosting_add-noneur': cb[6],
        'boosting_nonadd-noneur': cb[1],
        'snpboost': cb[8],
        'snpboost_noneur': cb[8],
        'snpboost_ex-chr6': cb[8],
        'snpnet': cb[5],
        'snpnet_nonadd': cb[5],  # cb[7],
        'snpnet_noneur': cb[5],
        'snpnet_ex-chr6': cb[5],
        'lassosum': cb[2],
        'lassosum_noneur': cb[2],
        'lassosum_ex-chr6': cb[2],
        'ldpred': cb[9],
        'ldpred_noneur': cb[9],
        'ldpred_ex-chr6': cb[9],
        'prscs': cb[0],
        'prscs_noneur': cb[0],
        'prscs_ex-chr6': cb[0],
        'sbayesr': cb[4],
        'sbayesr_noneur': cb[4],
        'sbayesr_ex-chr6': cb[4],
        'clump': cb[7],
        'clump_noneur': cb[7],
        'clump_ex-chr6': cb[7],
    }

    # methods = 'boosting_integrate,boosting_add,boosting_nonadd,snpnet,snpnet_nonadd,lassosum,ldpred,prscs,sbayesr,clump'.split(',')
    # mypalette = [
    #    cb[3], cb[1], cb[6],  # boosting red, orange, pink
    #    cb[8],  # yellow #lasso
    #    cb[7],
    #    cb[2],  # green # lassosum
    #    cb[9], cb[0],  # blue #ldpred, prs-cs
    #    cb[4],  # purple sbayesr
    #    cb[5],  # C+T
    # ]
    # assert len(mypalette) == len(methods)
    # paletted = dict(zip(methods, mypalette))
    return paletted


def accd():
    return {'auc': 'AUC',
            'pauc_fpr0.1': 'pAUC',
            'pauc_fpr0.05': 'pAUC',
            'auprc': 'AUPRC',
            'nagelkerke': 'Nagelkerke\'s $R^2$',
            'adj_mcfadden': 'Adjusted McFadden $R^2$'}


def statd():
    return {'jaccard_similarity': 'Jaccard Similarity',
            'chisq': 'Chi-square Test'}


# def methodd_capital(methods):
#    methodd = _methodd_capital()
#    check_keys_in_list(methodd, methods)
#    return methodd


def check_keys_in_list(d, li):
    for v in li:
        if v not in d:
            raise RuntimeError('Key not in dictionay: ', v, d)


def methodd(methods, kind=None):
    if kind is None:
        methodd = dict(zip(methods, methods))
        return methodd
    else:
        if kind == 'capital':
            methodd = _methodd_capital()
        elif kind == 'capital2':
            methodd = _methodd_capital2()
        elif kind == 'capital_roman':
            methodd = _methodd_capital_roman(methods,
                                             _create_xtick_labels_roman_number, _methodd_capital)
            # methodd = _methodd_capital_roman(methods)
        elif kind == 'capital2_roman':
            methodd = _methodd_capital_roman(methods,
                                             _create_xtick_labels_roman_number, _methodd_capital2)
            # methodd = _methodd_capital2_roman(methods)
        elif kind == 'capital_roman_serial':
            # number depends on order of methods
            methodd = _methodd_capital_roman(methods,
                                             _create_xtick_labels_roman_number_serial, _methodd_capital2)
            # methodd = _methodd_capital_roman_serial(methods)
        # elif kind == 'capital_roman_boosting':
        #    methodd = _methodd_capital_roman_boosting_unspecific(methods)
        # elif kind == 'capital_roman_boosting_snpnet':
        #    methodd = _methodd_capital_roman_boosting_snpnet_unspecific(methods)
        else:
            raise RuntimeError('Unknown kind', kind)
        check_keys_in_list(methodd, methods)
        return methodd


# every roman number is specified
def _methodd_capital_roman(methods, roman_fn, methodd_fn):
    roman = roman_fn(methods)
    # roman = _create_xtick_labels_roman_number(methods)
    methodd_cap = methodd_fn()
    # methodd_cap = _methodd_capital()
    methodd = {k: tick + '. ' + methodd_cap[k] for k, tick in zip(methods, roman)}
    return methodd


#  merged with _methodd_capital_roman()
# def _methodd_capital_roman_serial(methods):
#    roman = _create_xtick_labels_roman_number_serial(methods)
#    methodd_cap = _methodd_capital()
#    methodd = {k: tick + '. ' + methodd_cap[k] for k, tick in zip(methods, roman)}
#    return methodd


# unspecific: roman number depends on methods order
# def _methodd_capital_roman_boosting_unspecific(methods):
#    roman = _create_xtick_labels_roman_number_boosting_unspecific(methods)
#    methodd_cap = _methodd_capital()
#    methodd = {k: tick + '. ' + methodd_cap[k] for k, tick in zip(methods, roman)}
#    return methodd

# def _methodd_capital_roman_boosting_snpnet_unspecific(methods):
#    roman = _create_xtick_labels_roman_number_boosting_snpnet(methods)
#    methodd_cap = _methodd_capital()
#    methodd = {k: tick + '. ' + methodd_cap[k] for k, tick in zip(methods, roman)}
#    return methodd


def _roman_numbers():
    strs = ['i', 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii', 'ix', 'x', 'xi', 'xii', 'xiii']
    return strs


def _roman_number(i):
    strs = _roman_numbers()
    return strs[i]


def _create_xtick_labels_roman_number_serial(n):
    strs = _roman_numbers()
    strs = strs[:n]
    return strs


def create_xtick_labels(methods=None, kind=None, n=None):
    if kind is None:
        if n is not None:
            labels = list(range(n))
        elif methods is not None:
            labels = methods
        else:
            raise RuntimeError
    else:
        if kind == 'roman_number':
            labels = _create_xtick_labels_roman_number(methods)
        elif kind == 'roman_number_serial':
            # i, ii, etc.
            labels = _create_xtick_labels_roman_number_serial(n)
        # elif kind == 'roman_number_boosting':
        #    labels = _create_xtick_labels_roman_number_boosting_unspecific(methods)
        # elif kind == 'roman_number_boosting_snpnet':
        #    labels = _create_xtick_labels_roman_number_boosting_snpnet(methods)
        else:
            raise RuntimeError('Unknown kind', kind)
    return labels


def _create_xtick_labels_roman_number(methods):
    romand = _romand_method()
    check_keys_in_list(romand, methods)
    labels = [romand[x] for x in methods]
    return labels


# def _create_xtick_labels_roman_number_boosting_unspecific(methods):
#    '''
#    'boosting_integrate' -> 'i'
#    'boosting_add' -> 'a'
#    'boosting_nonadd' -> 'b'
#    'snpnet' -> 'ii'
#    '''
#    labels = []
#    ct = 0
#    ct_alphabet = 0
#    for method in methods:
#        if method != 'boosting_integrate' and method.startswith('boosting_'):
#            label = chr(ct_alphabet + ord('a'))
#            ct_alphabet += 1
#        else:
#            label = _roman_number(ct)
#            ct += 1
#        labels.append(label)
#    return labels


# def _create_xtick_labels_roman_number_boosting_snpnet(methods):
#    '''
#    'boosting_integrate' -> 'i'
#    'boosting_add' -> 'a'
#    'boosting_nonadd' -> 'b'
#    'snpnet' -> 'ii'
#    'snpnet-nonadd' -> 'iii'
#    'lassosum' -> 'iii'
#    '''
#    labels = []
#    ct = 0
#    ct_alphabet = 0
#    for method in methods:
#        if method != 'boosting_integrate' and method.startswith('boosting_'):
#            label = chr(ct_alphabet + ord('a'))
#            ct_alphabet += 1
#        elif method == 'snpnet_nonadd':
#            label = 'ii-a'
#        else:
#            label = _roman_number(ct)
#            ct += 1
#        labels.append(label)
#    return labels


def palette_hls_2_methods_strip_nsnv():
    # for marker palette
    # make some methods' color same as hls_2

    paletted = palette_hls_2_methods()

    # overwrite black
    methods_black = 'boosting_integrate,boosting_add,boosting_nonadd,snpboost,snpnet,snpnet_nonadd,lassosum,clump'.split(',')
    paletted |= {k: 'k' for k in methods_black}

    return paletted


def palette_diesase_status():
    # labels = ['Case', 'Control']
    labels = ['case', 'control']
    cb = sns.color_palette("colorblind", 10)
    colors = [cb[1], cb[0]]

    return dict(zip(labels, colors))


def phed(phes, kind=None):
    if kind is None:
        phed = dict(zip(phes, phes))
        return phed
    else:
        if kind == 'capital':
            phed = _phed_capital()
        elif kind == 'capital2':
            phed = _phed_capital2()
        elif kind == 'capital3':
            phed = _phed_capital3()
        else:
            raise RuntimeError('Unknown kind', kind)
        check_keys_in_list(phed, phes)
        return phed


def _phed_capital():
    phed = {
        't2d': 'Type 2 diabetes',
        'cad': 'Coronary artery disease',
        'bc': 'Breast cancer',
        'af': 'Atrial fibrillation',
        'ibd': 'Inflammatory bowel disease',
        'gout': 'Gout',
        'ra': 'Rheumatoid arthritis',
        'acd': 'All-cause dementia',
        'ad': 'Alzheimer\'s disease',
        'psr': 'Psoriasis',
        'cc': 'Colorectal cancer',
        'atm': 'Asthma',
        'her': 'Hernia',
        'ingher': 'Inguinal Hernia',
        'ami': 'Acute myocardial infarction',
        'mch': 'Major coronary heart disease event',
        'mi': 'Myocardial infarction',
        'hae': 'Hemorrhoids',
        'mbe': 'Myopia: Both eyes'
    }
    return phed


def _phed_capital2():
    phed = {
        't2d': 'Type 2 diabetes',
        'cad': 'Coronary artery disease',
        'bc': 'Breast cancer',
        'af': 'Atrial fibrillation',
        'ibd': 'Inflammatory bowel disease',
        'gout': 'Gout',
        'ra': 'Rheumatoid arthritis',
        'acd': 'All-cause dementia',
        'ad': 'Alzheimer\'s disease',
        'psr': 'Psoriasis',
        'cc': 'Colorectal cancer',
        'atm': 'Asthma',
        'her': 'Hernia',
        'ingher': 'Inguinal Hernia',
        'ami': 'Acute myocardial infarction',
        'mch': 'Major coronary\nheart disease event',
        'mi': 'Myocardial infarction',
        'hae': 'Hemorrhoids',
        'mbe': 'Myopia: Both eyes'
    }
    return phed


def _phed_capital3():
    phed = {
        't2d': 'Type 2\ndiabetes',
        'cad': 'Coronary artery\ndisease',
        'bc': 'Breast cancer',
        'af': 'Atrial\nfibrillation',
        'ibd': 'Inflammatory\nbowel disease',
        'gout': 'Gout',
        'ra': 'Rheumatoid\narthritis',
        'acd': 'All-cause\ndementia',
        'ad': 'Alzheimer\'s\ndisease',
        'psr': 'Psoriasis',
        'cc': 'Colorectal\ncancer',
        'atm': 'Asthma',
        'her': 'Hernia',
        'ingher': 'Inguinal Hernia',
        'ami': 'Acute myocardial\ninfarction',
        'mch': 'Major coronary\nheart disease event',
        'mi': 'Myocardial\ninfarction',
        'hae': 'Hemorrhoids',
        'mbe': 'Myopia:\nBoth eyes'
    }
    return phed


def phes_12():
    # TODO: automatically order by auc?
    return 'ra,psr,gout,ibd,atm,acd,ad,af,bc,cc,cad,t2d'.split(',')
    #return 'ra,atm,psr,gout,ibd,bc,acd,ad,cc,af,cad,t2d'.split(',')
    # return 't2d,cad,bc,af,ibd,ra,gout,ad,acd,cc,atm,psr'.split(',')


def statusd():
    # return {0: 'Control', 1: 'Case'}
    return {'control': 'Control', 'case': 'Case'}


def genetic_models():
    # order matters
    return ['overdom', 'dom', 'add', 'rec', 'overrec']


def genetic_modeld():
    return {'overrec': 'Overrecessive', 'rec': 'Recessive', 'add': 'Additive', 'dom': 'Dominant', 'overdom': 'Overdominant'}


def genetic_modeld2():
    return {'overrec': 'Overrec', 'rec': 'Rec', 'add': 'Add', 'dom': 'Dom', 'overdom': 'Overdom'}


def palette_hls_2_genetic_model():

    cb = sns.color_palette("colorblind", 10)
    paletted = {'overdom': cb[6],
                'dom': cb[4],
                'add': cb[2],
                'rec': cb[0],
                'overrec': cb[9]
                }
    return paletted
