$include "data.gms"
$onEmpty
$onInline

Parameters 
    Dem_ext(k, t)
    Ret_ext(k, t)
    RemanFr_ext(k, t);

Binary Variables
    WI(i)                                           1 if a distribution center is built at site j
    
Positive Variables
    XI(i)                                           Capacity of distribution center built at site j
    PRI(i, t)                                        New product produced at site i in time period t
    FIJ(i, j, t)                                    Flow from site i to j in time period t
    FJK(j, k, t)                                    Flow from site j to customer k in time period t
    FKC(k, c, t)                                    Flow from customer k to site c in time period t
    FCI(c, i, t)                                    Flow from site c to site i in time period t
    FCD(c, d, t)                                    Flow from site c to disposal center d in time period t
    
    BSI(i, t)                                       Basestock level in distribution center at site j;

Integer Variables
    ZJ(j, t)                                        Number of production modules constructed at site i in time period t
    ZC(c, t)                                        Number of collection modules constructed at site c in time period t
    YJ(j, j_, t)                                    Number of modules transferred from i to i_ in time period t
    YC(c, c_, t)                                    Number of modules transferred from c to c_ in time period t
    VJ(j, t)                                        Number of production modules present at site i in time period t
    VC(c, t)                                        Number of collection modules present at site c in time period t
    
Variable
    obj_MCLSCN;
    
Equations
    obj_cons
    prod_construct_cons1(i)
    prod_construct_cons2(i)
    demand_cons(k, t)
    return_cons(k, t)
    
    prod_capacity_cons(i, t)
    dist_capacity_cons(j, t)
    coll_capacity_cons(c, t)
    
    dist_module_limit_cons(j, t)
    coll_module_limit_cons(c, t)
    dist_module_count_cons(j, t)
    coll_module_count_cons(c, t)

    flow_balance_i_cons(i, t)
    flow_balance_j_cons(j, t)
    flow_balance_c_cons(c, t)
    reman_fraction_cons(c, t)
    inventory_lim_cons(i, t);
    
obj_cons ..                                     obj_MCLSCN =g= sum(i, BaseCostI(i) * WI(i) + CapCostI(i) * XI(i) + sum(t, InvCostI(i) * BSI(i, t)))
                                                                + sum(j, sum(t, ModCostJ * ZJ(j, t) + sum(j_ $ (not sameAs(j, j_)), RelocCostJ(j, j_) * YJ(j, j_, t))))
                                                                + sum(c, sum(t, ModCostC * ZC(c, t) + sum(c_ $ (not sameAs(c, c_)), RelocCostC(c, c_) * YC(c, c_, t))))
                                                                + sum((i, t), ManCost(i) * PRI(i, t))
                                                                + sum((c, i, t), RemanCost(i) * FCI(c, i, t))
                                                                + sum((k, c, t), ColCost(c) * FKC(k, c, t))
                                                                + sum((c, d, t), DispCost(d) * FCD(c, d, t))
                                                                + sum((i, j, t), TrCostIJ(i, j) * FIJ(i, j, t))
                                                                + sum((j, k, t), TrCostJK(j, k) * FJK(j, k, t))
                                                                + sum((k, c, t), TrCostKC(k, c) * FKC(k, c, t))
                                                                + sum((c, d, t), TrCostCD(c, d) * FCD(c, d, t))
                                                                + sum((c, i, t), TrCostCI(c, i) * FCI(c, i, t))
                                                                + sum(t, sum(j, VJ(j, t) * 2000) + sum(c, VC(c, t) * 1000));
prod_construct_cons1(i)  ..                     XI(i) =l= CapMaxI * WI(i);
prod_construct_cons2(i)  ..                     XI(i) =g= CapMinI * WI(i);
demand_cons(k, t) ..                            sum(j, FJK(j, k, t)) =g= Dem_ext(k, t);
return_cons(k, t) ..                            sum(c, FKC(k, c, t)) =e= Ret_ext(k, t);
prod_capacity_cons(i, t) ..                     BSI(i, t - 1) + PRI(i, t) + sum(c, FCI(c, i, t)) =l= XI(i); 
dist_capacity_cons(j, t) ..                     sum(k, FJK(j, k, t)) =l= VJ(j, t) * ModSizeJ;
coll_capacity_cons(c, t) ..                     sum(k, FKC(k, c, t)) =l= VC(c, t) * ModSizeC;
dist_module_limit_cons(j, t) ..                 VJ(j, t) =l= MaxModNumJ;
coll_module_limit_cons(c, t) ..                 VC(c, t) =l= MaxModNumC;
dist_module_count_cons(j, t) ..                 VJ(j, t) =e= VJ(j, t - 1) + ZJ(j, t) - sum(j_ $ (not sameAs(j, j_)), YJ(j, j_, t)) + sum(j_ $ (not sameAs(j, j_)), YJ(j_, j, t));
coll_module_count_cons(c, t) ..                 VC(c, t) =e= VC(c, t - 1) + ZC(c, t) - sum(c_ $ (not sameAs(c, c_)), YC(c, c_, t)) + sum(c_ $ (not sameAs(c, c_)), YC(c_, c, t));
flow_balance_i_cons(i, t) ..                    BSI(i, t - 1) + PRI(i, t) + sum(c, FCI(c, i, t)) =e= sum(j, FIJ(i, j, t)) + BSI(i, t);
flow_balance_j_cons(j, t) ..                    sum(i, FIJ(i, j, t)) =e= sum(k, FJK(j, k, t));
flow_balance_c_cons(c, t) ..                    sum(k, FKC(k, c, t)) =e= sum(i, FCI(c, i, t)) + sum(d, FCD(c, d, t));
reman_fraction_cons(c, t) ..                    sum(i, FCI(c, i, t)) =l= sum(k, RemanFr_ext(k, t) * FKC(k, c, t));
inventory_lim_cons(i, t) ..                     BSI(i, t) =l= XI(i);

ZJ.up(j, t) = MaxModNumJ;
ZC.up(c, t) = MaxModNumC;



Model mod_MCLSCN_det /
    obj_cons
    prod_construct_cons1
    prod_construct_cons2
    demand_cons
    return_cons
    prod_capacity_cons
    dist_capacity_cons
    coll_capacity_cons
    dist_module_limit_cons
    coll_module_limit_cons
    dist_module_count_cons
    coll_module_count_cons
    flow_balance_i_cons
    flow_balance_j_cons
    flow_balance_c_cons
    reman_fraction_cons
    inventory_lim_cons
/;
Option mip = gurobi;
Option threads = -1;
Option optCR = 0;

Set
    smp                 / 1 * 200 /;

Parameters
    Transport_cost
    Process_cost
    Construct_cost
    ObjWReloc(smp)
    ObjWOReloc(smp)
    YJ_rec(j, j_, t)
    YC_rec(c, c_, t);
    
$onText
Transport_cost = sum((i, j, t), TrCostIJ(i, j) * FIJ.l(i, j, t))
                + sum((j, k, t), TrCostJK(j, k) * FJK.l(j, k, t))
                + sum((k, c, t), TrCostKC(k, c) * FKC.l(k, c, t))
                + sum((c, d, t), TrCostCD(c, d) * FCD.l(c, d, t))
                + sum((c, i, t), TrCostCI(c, i) * FCI.l(c, i, t));
Process_cost = sum((i, t), ManCost(i) * PRI.l(i, t))
                + sum((c, i, t), RemanCost(i) * FCI.l(c, i, t))
                + sum((k, c, t), ColCost(c) * FKC.l(k, c, t))
                + sum((c, d, t), DispCost(d) * FCD.l(c, d, t));
Construct_cost = sum(i, BaseCostI(i) * WI.l(i) + CapCostI(i) * XI.l(i) + sum(t, InvCostI(i) * BSI.l(i, t)))
                + sum(j, sum(t, ModCostJ * ZJ.l(j, t) + sum(j_ $ (not sameAs(j, j_)), RelocCostJ(j, j_) * YJ.l(j, j_, t))))
                + sum(c, sum(t, ModCostC * ZC.l(c, t) + sum(c_ $ (not sameAs(c, c_)), RelocCostC(c, c_) * YC.l(c, c_, t))));
$offText

execseed = 1 + gmillisec(jnow);
loop(smp, 
    loop(t, 
        Dem_ext(k, t) = sum(s $ (ord(s) = uniformInt(1, card(s))), Dem(k, t, s));
        Ret_ext(k, t) = sum(s $ (ord(s) = uniformInt(1, card(s))), Ret(k, t, s));
        RemanFr_ext(k, t) = sum(s $ (ord(s) = uniformInt(1, card(s))), RemanFr(k, t, s));
    );
    Solve mod_MCLSCN_det using mip minimizing obj_MCLSCN;
    ObjWReloc(smp) = obj_MCLSCN.l;
    YJ_rec(j, j_, t) = YJ.l(j, j_, t);
    YC_rec(c, c_, t) = YC.l(c, c_, t);

    YJ.up(j, j_, t) = 0;
    YC.up(c, c_, t) = 0;
    Solve mod_MCLSCN_det using mip minimizing obj_MCLSCN;
    ObjWOReloc(smp) = obj_MCLSCN.l;
);
execute_unload "MCLSCN_det_result.gdx";
