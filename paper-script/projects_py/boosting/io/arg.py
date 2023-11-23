
from collections import defaultdict


def arg_method_mid_path(method_model_mid_path_arg):
    # method-> mid_path
    # tuple not list for @lru.cache()
    # ex. (kind, ), (kind, model), (kind, simu_type, model)
    method_mid_path = defaultdict(lambda: ('default',))
    # TODO: if method_model_mid_path_arg is None: return
    for v in method_model_mid_path_arg:

        if len(v) < 2:
            raise RuntimeError('Unknown arg len in --method-model: ', v)

        method, mid_path = v[0], tuple(v[1:])
        d = {method: mid_path}
        method_mid_path |= d

        # if len(v) == 2:
        #    # [method, kind]
        #    method, kind = v[0], v[1]
        #    d = {method: (kind, None)}
        #    method_modeld |= d
        #    # method_modeld[method]=(kind,None)
        # elif len(v) == 3:
        #    # [method, kind, model]
        #    method, kind, model = v[0], v[1], v[2]
        #    d = {method: (kind, model)}
        #    method_modeld |= d
        # else:
        # raise RuntimeError('Unknown arg len in --method-model: ', v)

    # method_modeld |= dict(args.method_model)
    # method_modeld |= {
    #    'boosting_nonadd': args.model_boost,
    #    'boosting_add': args.model_boost_add,
    #    'boosting_loss': args.model_boost_loss
    # }
    # logger.debug('method_modeld {}'.format(method_modeld))
    return method_mid_path


def arg_method_type(method_model_type_arg):
    # method-> type (e.g. noneur, ex-chr6)

    methodd_type = defaultdict(lambda: None)
    if method_model_type_arg is None:
        return methodd_type

    if len(method_model_type_arg) < 2:
        raise RuntimeError('Unknown arg len in --method-model: ', method_model_type_arg)

    mtype, methods = method_model_type_arg[0], list(method_model_type_arg[1:])
    for method in methods:
        d = {method: mtype}
        methodd_type |= d

    return methodd_type
