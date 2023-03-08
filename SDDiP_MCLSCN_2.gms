$include "data_2.gms"
$onEmpty
$onInline

Set
    bin                                             Set for binary dissociation      / 1 * 3 /
    m                                               Set for samples per iteration    / 1 * %mcount% /
    iter                                            Set for iterations               / 1 * 200 /;
    
Alias (bin, bin_);
Alias (m, m_);

Parameters
*** Powers of 2 for binary decomposition for integer parameters
    binDigits(bin)
    
*** Random sample for each stage (elements take 1..s)
    scenarios(m, t)
    
*** Punishment for demand slack
    DemPunish                               / 8000 /

*** Values from previous stage
    VJ_prev(j, bin)                         / /
    VC_prev(c, bin)                         / /

*** Collect solutions from forward problems
    YJ_fsol(j, j_, t, m, iter)              / /
    YC_fsol(c, c_, t, m, iter)              / /
    VJ_fsol(j, bin, t, m, iter)             / /
    VC_fsol(c, bin, t, m, iter)             / /
    theta_fsol(t, m, iter)                  / /
    
*** Time stage & scenario-related parameters
    Dem_ext(k)
    Ret_ext(k)
    RemanFr_ext(k)
    ModCostJ_ext(j)
    ModCostC_ext(c)
    
    sigma_b(m_, s_, iter)                   / /
    pi_b_VJ(j, bin, m_, s_, iter)           / /
    pi_b_VC(c, bin, m_, s_, iter)           / /
    PRI_cutparam(i, m_, s_, iter)           / /
    VJ_cutparam(j, bin, m_, iter)           / /
    VC_cutparam(c, bin, m_, iter)           / /
    sigma_i(m_, s_, iter)                   / /
    
*** Cut parameters
    sigma_b_stages(t, m_, s_, iter)             / /
    pi_b_VJ_stages(j, bin, t, m_, s_, iter)     / /
    pi_b_VC_stages(c, bin, t, m_, s_, iter)     / /
    VJ_cutparam_stages(j, bin, t, m_, iter)     / /
    VC_cutparam_stages(c, bin, t, m_, iter)     / /
    sigma_i_stages(t, m_, s_, iter)             / /
    
*** Objective collection parameters
    UB(iter)                                / /
    LB(iter)                                / /
    um(m, t)
    mean
    mean_record(iter)
    stdev
    stdev_record(iter);    

Scalar runTime;
Scalar itc, t_idx, objOP, subgr_idx;

runTime = gsecond(jnow);
binDigits(bin) = power(2, ord(bin) - 1);

Positive Variables
    PRI_fw(i)                                       New product produced at site i in time period t
    FIJ_fw(i, j)                                    Flow from site i to j in time period t
    FJK_fw(j, k)                                    Flow from site j to customer k in time period t
    FKC_fw(k, c)                                    Flow from customer k to site c in time period t
    FCI_fw(c, i)                                    Flow from site c to site i in time period t
    FCD_fw(c, d)                                    Flow from site c to disposal center d in time period t
    SLKDem_fw(k)
    SLKRet_fw(k)
    
    VJ_dup(j, bin)
    VC_dup(c, bin)
    theta;

Integer Variables
    ZJ_fw(j)                                        Number of production modules constructed at site i in time period t
    ZC_fw(c)                                        Number of collection modules constructed at site c in time period t
    YJ_fw(j, j_)                                    Number of modules transferred from i to i_ in time period t
    YC_fw(c, c_)                                    Number of modules transferred from c to c_ in time period t
    
Binary Variables
    VJ_fw(j, bin)                                   Number of production modules present at site i in time period t
    VC_fw(c, bin)                                   Number of collection modules present at site c in time period t
    
Variable
    obj_fw;
    
Equations
    obj_cons
    demand_cons(k)
    return_cons(k)
    
    prod_capacity_cons(i)
    dist_capacity_cons(j)
    coll_capacity_cons(c)
    
    dist_module_limit_cons(j)
    coll_module_limit_cons(c)
    dist_zlimit_cons(j)
    coll_zlimit_cons(c)
    dist_module_count_cons(j)
    coll_module_count_cons(c)

    flow_balance_i_cons(i)
    flow_balance_j_cons(j)
    flow_balance_c_cons(c)
    reman_fraction_cons(c)
    
    VJ_link_cons(j, bin)
    VC_link_cons(c, bin)
    
    sbenders_cut_cons(m_, iter)
    integer_cut_cons(m_, iter);
    
obj_cons ..                                         obj_fw =g= sum(j, ModCostJ_ext(j) * ZJ_fw(j) + sum(j_ $ (not sameAs(j, j_)), RelocCostJ(j, j_) * YJ_fw(j, j_)))
                                                                    + sum(c, ModCostC_ext(c) * ZC_fw(c) + sum(c_ $ (not sameAs(c, c_)), RelocCostC(c, c_) * YC_fw(c, c_)))
                                                                    + sum(i, ManCost(i) * PRI_fw(i))
                                                                    + sum((c, i), RemanCost(i) * FCI_fw(c, i))
                                                                    + sum((k, c), ColCost(c) * FKC_fw(k, c))
                                                                    + sum((c, d), DispCost(d) * FCD_fw(c, d))
                                                                    + sum((i, j), TrCostIJ(i, j) * FIJ_fw(i, j))
                                                                    + sum((j, k), TrCostJK(j, k) * FJK_fw(j, k))
                                                                    + sum((k, c), TrCostKC(k, c) * FKC_fw(k, c))
                                                                    + sum((c, d), TrCostCD(c, d) * FCD_fw(c, d))
                                                                    + sum((c, i), TrCostCI(c, i) * FCI_fw(c, i))
                                                                    + theta;
demand_cons(k) ..                                   sum(j, FJK_fw(j, k)) =g= Dem_ext(k);
return_cons(k) ..                                   sum(c, FKC_fw(k, c)) =e= Ret_ext(k);
prod_capacity_cons(i) ..                            PRI_fw(i) + sum(c, FCI_fw(c, i)) =l= CapI(i); 
dist_capacity_cons(j) ..                            sum(k, FJK_fw(j, k)) =l= sum(bin, binDigits(bin) * VJ_fw(j, bin)) * ModSizeJ;
coll_capacity_cons(c) ..                            sum(k, FKC_fw(k, c)) =l= sum(bin, binDigits(bin) * VC_fw(c, bin)) * ModSizeC;
dist_module_limit_cons(j) ..                        sum(bin, binDigits(bin) * VJ_fw(j, bin)) =l= MaxModNumJ;
coll_module_limit_cons(c) ..                        sum(bin, binDigits(bin) * VC_fw(c, bin)) =l= MaxModNumC;
dist_zlimit_cons(j) ..                              ZJ_fw(j) =l= MaxModNumJ;
coll_zlimit_cons(c) ..                              ZC_fw(c) =l= MaxModNumC;
dist_module_count_cons(j) ..                        sum(bin, binDigits(bin) * VJ_fw(j, bin)) =e= sum(bin, binDigits(bin) * VJ_dup(j, bin)) + ZJ_fw(j) - sum(j_ $ (not sameAs(j, j_)), YJ_fw(j, j_)) + sum(j_ $ (not sameAs(j, j_)), YJ_fw(j_, j));
coll_module_count_cons(c) ..                        sum(bin, binDigits(bin) * VC_fw(c, bin)) =e= sum(bin, binDigits(bin) * VC_dup(c, bin)) + ZC_fw(c) - sum(c_ $ (not sameAs(c, c_)), YC_fw(c, c_)) + sum(c_ $ (not sameAs(c, c_)), YC_fw(c_, c));
flow_balance_i_cons(i) ..                           PRI_fw(i) + sum(c, FCI_fw(c, i)) =e= sum(j, FIJ_fw(i, j));
flow_balance_j_cons(j) ..                           sum(i, FIJ_fw(i, j)) =e= sum(k, FJK_fw(j, k));
flow_balance_c_cons(c) ..                           sum(k, FKC_fw(k, c)) =e= sum(i, FCI_fw(c, i)) + sum(d, FCD_fw(c, d));
reman_fraction_cons(c) ..                           sum(i, FCI_fw(c, i)) =l= sum(k, RemanFr_ext(k) * FKC_fw(k, c));
VJ_link_cons(j, bin) ..                             VJ_dup(j, bin) =e= VJ_prev(j, bin);
VC_link_cons(c, bin) ..                             VC_dup(c, bin) =e= VC_prev(c, bin);
sbenders_cut_cons(m_, iter) ..
                                                    theta =g= sum(s_, sigma_b(m_, s_, iter) + sum((j, bin), pi_b_VJ(j, bin, m_, s_, iter) * VJ_fw(j, bin)) + sum((c, bin), pi_b_VC(c, bin, m_, s_, iter) * VC_fw(c, bin))) / card(s_);
integer_cut_cons(m_, iter) ..
                                                    theta =g= (sum(s_, sigma_i(m_, s_, iter) * (sum((j, bin), (VJ_cutparam(j, bin, m_, iter) - 1) * VJ_fw(j, bin) + VJ_cutparam(j, bin, m_, iter) * (VJ_fw(j, bin) - 1)) + sum((c, bin), (VC_cutparam(c, bin, m_, iter) - 1) * VC_fw(c, bin) + VC_cutparam(c, bin, m_, iter) * (VC_fw(c, bin) - 1)))) + sum(s_, sigma_i(m_, s_, iter))) / card(s_);

*YJ_fw.fx(j, j_) = 0;
*YC_fw.fx(c, c_) = 0;

Model mod_fw /
    obj_cons
    demand_cons
    return_cons
    prod_capacity_cons
    dist_capacity_cons
    coll_capacity_cons
    dist_module_limit_cons
    coll_module_limit_cons
    dist_zlimit_cons
    coll_zlimit_cons
    dist_module_count_cons
    coll_module_count_cons
    flow_balance_i_cons
    flow_balance_j_cons
    flow_balance_c_cons
    reman_fraction_cons
    VJ_link_cons
    VC_link_cons
    sbenders_cut_cons
    integer_cut_cons
/;

Parameter
    pi_lag_VJ(j, bin)
    pi_lag_VC(c, bin)
    pi_lag_VJ_all(j, bin, m, s)
    pi_lag_VC_all(c, bin, m, s)
    
    LambdaMW                                        / 0.8 /
    VJ_core(j, bin)                                 / /
    VC_core(c, bin)                                 / /
    VJ_core_all(j, bin, t, m)                       / /
    VC_core_all(c, bin, t, m)                       / /;

Variable
    obj_lag;

Equation
    obj_lagrangean_cons
    VJ_link_MW_cons(j, bin)
    VC_link_MW_cons(c, bin);
    
obj_lagrangean_cons ..                              obj_lag =e= obj_fw - sum((j, bin), pi_lag_VJ(j, bin) * VJ_dup(j, bin)) - sum((c, bin), pi_lag_VC(c, bin) * VC_dup(c, bin));
VJ_link_MW_cons(j, bin) ..                          VJ_dup(j, bin) =e= VJ_core(j, bin);
VC_link_MW_cons(c, bin) ..                          VC_dup(c, bin) =e= VC_core(c, bin);

Model mod_imw /
    mod_fw - VJ_link_cons - VC_link_cons + VJ_link_MW_cons + VC_link_MW_cons
/;

Model mod_lag /
    mod_fw - VJ_link_cons - VC_link_cons + obj_lagrangean_cons
/;

VJ_dup.up(j, bin) = 1;
VC_dup.up(c, bin) = 1;

execseed = 1 + gmillisec(jnow);
Option optCR = 0;
Option mip = gurobi;
Option rmip = gurobi;
Option threads = -1;

Scalar algGap, LBimprove;
Scalar endFlag / 0 /; 

for (itc = 1 to card(iter),
    scenarios(m, t) = uniformInt(1, card(s));
*** Forward step
    loop (t,
        ModCostJ_ext(j) = ModCostJ(j, t);
        ModCostC_ext(c) = ModCostC(c, t);
            
        sigma_b(m_, s_, iter) = sigma_b_stages(t, m_, s_, iter);
        pi_b_VJ(j, bin, m_, s_, iter) = pi_b_VJ_stages(j, bin, t, m_, s_, iter);
        pi_b_VC(c, bin, m_, s_, iter) = pi_b_VC_stages(c, bin, t, m_, s_, iter);
        VJ_cutparam(j, bin, m_, iter) = VJ_fsol(j, bin, t, m_, iter);
        VC_cutparam(c, bin, m_, iter) = VC_fsol(c, bin, t, m_, iter);
        sigma_i(m_, s_, iter) = sigma_i_stages(t, m_, s_, iter);
        loop(m,
            Dem_ext(k) = sum(s $ (ord(s) = round(scenarios(m, t))), Dem(k, t, s));
            Ret_ext(k) = sum(s $ (ord(s) = round(scenarios(m, t))), Ret(k, t, s));
            RemanFr_ext(k) = sum(s $ (ord(s) = round(scenarios(m, t))), RemanFr(k, t, s));
            
            VJ_prev(j, bin) = sum(iter $ (ord(iter) = itc), VJ_fsol(j, bin, t - 1, m, iter) $ (ord(t) > 1));
            VC_prev(c, bin) = sum(iter $ (ord(iter) = itc), VC_fsol(c, bin, t - 1, m, iter) $ (ord(t) > 1));
            
            Solve mod_fw using mip minimizing obj_fw;
            
            VJ_fsol(j, bin, t, m, iter) $ (ord(iter) = itc) = round(VJ_fw.l(j, bin));
            VC_fsol(c, bin, t, m, iter) $ (ord(iter) = itc) = round(VC_fw.l(c, bin));
            YJ_fsol(j, j_, t, m, iter) $ (ord(iter) = itc) = round(YJ_fw.l(j, j_));
            YC_fsol(c, c_, t, m, iter) $ (ord(iter) = itc) = round(YC_fw.l(c, c_));
            theta_fsol(t, m, iter) $ (ord(iter) = itc) = theta.l;
            um(m, t) = obj_fw.l - theta.l;
            
*            if (ord(m) = 1,
*                execute_unload "forw_intervdata.gdx";
*            );
        );
    );
    
*** Statistical upper bound update
    mean = sum((m, t), um(m, t)) / card(m);
    mean_record(iter) $ (ord(iter) = itc) = mean;
    stdev = sqrt(sum(m, sqr(sum(t, um(m, t)) - mean)) / (card(m) - 1));
    stdev_record(iter) $ (ord(iter) = itc) = stdev;
    UB(iter) $ (ord(iter) = itc) = mean + 1.96 * stdev / sqrt(card(m));
    
*** Update core points
    if (itc = 1,
        VJ_core_all(j, bin, t, m) = sum(iter $ (ord(iter) = itc), VJ_fsol(j, bin, t - 1, m, iter));
        VC_core_all(c, bin, t, m) = sum(iter $ (ord(iter) = itc), VC_fsol(c, bin, t - 1, m, iter));
    else
        VJ_core_all(j, bin, t, m) = LambdaMW * sum(iter $ (ord(iter) = itc), VJ_fsol(j, bin, t - 1, m, iter)) + (1 - LambdaMW) * VJ_core_all(j, bin, t, m);
        VC_core_all(c, bin, t, m) = LambdaMW * sum(iter $ (ord(iter) = itc), VC_fsol(c, bin, t - 1, m, iter)) + (1 - LambdaMW) * VC_core_all(c, bin, t, m);
    );

*** Backward step
    for (t_idx = card(t) downto 2,
        ModCostJ_ext(j) = sum(t $ (ord(t) = 1), ModCostJ(j, t));
        ModCostC_ext(c) = sum(t $ (ord(t) = 1), ModCostC(c, t));
        sigma_b(m_, s_, iter) = sum(t $ (ord(t) = t_idx), sigma_b_stages(t, m_, s_, iter));
        pi_b_VJ(j, bin, m_, s_, iter) = sum(t $ (ord(t) = t_idx), pi_b_VJ_stages(j, bin, t, m_, s_, iter));
        pi_b_VC(c, bin, m_, s_, iter) = sum(t $ (ord(t) = t_idx), pi_b_VC_stages(c, bin, t, m_, s_, iter));
        VJ_cutparam(j, bin, m_, iter) = sum(t $ (ord(t) = t_idx), VJ_fsol(j, bin, t, m_, iter));
        VC_cutparam(c, bin, m_, iter) = sum(t $ (ord(t) = t_idx), VC_fsol(c, bin, t, m_, iter));
        sigma_i(m_, s_, iter) = sum(t $ (ord(t) = t_idx), sigma_i_stages(t, m_, s_, iter));
        
****** LP Relaxation step of the Strengthened Benders cut
*$onText
        loop (m,
            VJ_prev(j, bin) = sum((t, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc), VJ_fsol(j, bin, t, m, iter));
            VC_prev(c, bin) = sum((t, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc), VC_fsol(c, bin, t, m, iter));
            VJ_core(j, bin) = sum(t $ (ord(t) = t_idx), VJ_core_all(j, bin, t, m));
            VC_core(c, bin) = sum(t $ (ord(t) = t_idx), VC_core_all(c, bin, t, m));
            
            loop(s,
                Dem_ext(k) = sum(t $ (ord(t) = t_idx), Dem(k, t, s));
                Ret_ext(k) = sum(t $ (ord(t) = t_idx), Ret(k, t, s));
                RemanFr_ext(k) = sum(t $ (ord(t) = t_idx), RemanFr(k, t, s));
                
                /* Strengthened Benders cut */
*$onText
                Solve mod_imw using rmip minimizing obj_fw;
                pi_b_VJ_stages(j, bin, t, m, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc) = VJ_link_MW_cons.m(j, bin);
                pi_b_VC_stages(c, bin, t, m, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc) = VC_link_MW_cons.m(c, bin);
                pi_lag_VJ(j, bin) = VJ_link_MW_cons.m(j, bin);
                pi_lag_VC(c, bin) = VC_link_MW_cons.m(c, bin);
*$offText
$onText            
                Solve mod_fw using rmip minimizing obj_fw;
                pi_b_VJ_stages(j, bin, t, m, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc) = VJ_link_cons.m(j, bin);
                pi_b_VC_stages(c, bin, t, m, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc) = VC_link_cons.m(c, bin);
                pi_lag_VJ(j, bin) = VJ_link_cons.m(j, bin);
                pi_lag_VC(c, bin) = VC_link_cons.m(c, bin);
$offText
                Solve mod_lag using mip minimizing obj_lag;
                sigma_b_stages(t, m, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc) = obj_lag.l;
                
                /* Integer cut */
                Solve mod_fw using mip minimizing obj_fw;
                sigma_i_stages(t, m, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc) = obj_fw.l;
            )
        );

****** Lagrangean cut
$onText
        mod_lag.solveLink = 2;
        loop((m, s),
            Dem_ext(k) = sum(t $ (ord(t) = t_idx), Dem(k, t, s));
            Ret_ext(k) = sum(t $ (ord(t) = t_idx), Ret(k, t, s));
            RemanFr_ext(k) = sum(t $ (ord(t) = t_idx), Ret(k, t, s));
            VJ_prev(j, bin) = sum((t, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc), VJ_fsol(j, bin, t, m, iter));
            VC_prev(c, bin) = sum((t, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc), VC_fsol(c, bin, t, m, iter));
            objOP = sum((t, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itc), sigma_i_stages(t, m, s, iter));
            
            for (subgr_idx = 1 to 5,
                Solve mod_lag using mip minimizing obj_lag;
                pi_lag_VJ(j, bin) = pi_lag_VJ(j, bin) - (objOP - obj_lag.l) / sqr(sum((j_, bin_), VJ_dup.l(j_, bin_)) + sum((c_, bin_), VC_dup.l(c_, bin_))) * (VJ_dup.l(j, bin) - VJ_prev(j, bin));
                pi_lag_VC(c, bin) = pi_lag_VC(c, bin) - (objOP - obj_lag.l) / sqr(sum((j_, bin_), VJ_dup.l(j_, bin_)) + sum((c_, bin_), VC_dup.l(c_, bin_))) * (VC_dup.l(c, bin) - VC_prev(c, bin));
            );
*            pi_l_stages(j, m, bin, t, k, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itercount) = pi_lag(j, m, bin);
*            sigma_l_stages(t, k, s, iter) $ (ord(t) = t_idx - 1 and ord(iter) = itercount) = obj_lag.l;
        );
$offText      
    );

*** Solve stage 1 problem to get LB
    VJ_prev(j, bin) = 0;
    VC_prev(c, bin) = 0;
    ModCostJ_ext(j) = sum(t $ (ord(t) = 1), ModCostJ(j, t));
    ModCostC_ext(c) = sum(t $ (ord(t) = 1), ModCostC(c, t));
    loop (s,
        Dem_ext(k) = sum(t $ (ord(t) = 1), Dem(k, t, s));
        Ret_ext(k) = sum(t $ (ord(t) = 1), Ret(k, t, s));
        RemanFr_ext(k) = sum(t $ (ord(t) = 1), RemanFr(k, t, s));
        
        sigma_b(m_, s_, iter) = sum(t $ (ord(t) = 1), sigma_b_stages(t, m_, s_, iter));
        pi_b_VJ(j, bin, m_, s_, iter) = sum(t $ (ord(t) = 1), pi_b_VJ_stages(j, bin, t, m_, s_, iter));
        pi_b_VC(c, bin, m_, s_, iter) = sum(t $ (ord(t) = 1), pi_b_VC_stages(c, bin, t, m_, s_, iter));
        VJ_cutparam(j, bin, m_, iter) = sum(t $ (ord(t) = 1), VJ_fsol(j, bin, t, m_, iter));
        VC_cutparam(c, bin, m_, iter) = sum(t $ (ord(t) = 1), VC_fsol(c, bin, t, m_, iter));
        sigma_i(m_, s_, iter) = sum(t $ (ord(t) = 1), sigma_i_stages(t, m_, s_, iter));

        Solve mod_fw using mip minimizing obj_fw;
        LB(iter) $ (ord(iter) = itc) = LB(iter) + obj_fw.l;
    );
    LB(iter) $ (ord(iter) = itc) = LB(iter) / card(s);
    
    execute_unload "SDDiP_MCLSCN_result.gdx";
    
    algGap = sum(iter $ (ord(iter) = itc), UB(iter) - LB(iter)) / sum(iter $ (ord(iter) = itc), UB(iter));
    if (itc > 10,
        LBimprove = sum(iter $ (ord(iter) = itc), LB(iter)) - sum(iter $ (ord(iter) = itc - 10), LB(iter));
        LBimprove = LBimprove / sum(iter $ (ord(iter) = itc), LB(iter));
    );
    if (algGap < 0.01,
        endFlag = 1;
    );
    if (((itc > 10) and (LBimprove < 1e-3)) or (endFlag > 0),
        break;
    )
);

runTime = timeExec;
execute_unload "SDDiP_MCLSCN_result.gdx";
