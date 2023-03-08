$onEmpty
Sets
    i                   Set for production sites
    j                   Set for distribution centers j
    k                   Set for customers
    c                   Set for collection and recycling sites
    d                   Set for disposal centers
    t                   Set for time stages
    s                   Set for scenarios per stage;
    
Alias (i, i_, i__);
Alias (j, j_, j__);
Alias (k, k_, k__);
Alias (c, c_, c__);
Alias (t, t_, t__);
Alias (s, s_, s__);

Parameters
    Dem(k, t, s)                                    Demand for customers
    Ret(k, t, s)                                    Returned product from customers
    RemanFr(k, t, s)                                Fraction that can be remanufacutured
    
    TrCostIJ(i, j)                                  Transportation cost from production site i to distribution site j
    TrCostJK(j, k)                                  Transportation cost from distribution site j to customer site k                                
    TrCostKC(k, c)                                  Transportation cost from production site i to distribution cen
    TrCostCD(c, d)                                  Transportation cost from production site i to distribution cen
    TrCostCI(c, i)                                  Transportation cost from production site i to distribution cen
    
    ModCostJ(j, t)                                  Cost for one distribution module
    ModCostC(c, t)                                  Cost for one collection & recycle module
    MaxModNumJ                                      Max number of modules installed at a distribution site
    MaxModNumC                                      Max number of modules installed at a collection site
    ModSizeJ                                        Size for per distribution module
    ModSizeC                                        Size for per collection & recycle module
    RelocCostJ(j, j_)                               Distribution module relocation cost from site j to j_
    RelocCostC(c, c_)                               Collection module relocation cost from site c to c_
    
    CapI(i)                                         Max capacity for production site i
    
    ManCost(i)                                      Unit manufacturing cost at production site i
    RemanCost(i)                                    Unit remanufacturing cost at production site i
    ColCost(c)                                      Unit collection & recycling cost at collection center c
    DispCost(d)                                     Unit disposal cost at disposal center d;
    
$gdxIn "data_2.gdx"
$load i, j, k, c, d, t, s
$load Dem, Ret, RemanFr, TrCostIJ, TrCostJK, TrCostKC, TrCostCD, TrCostCI
$load ModCostJ, ModCostC, MaxModNumJ, MaxModNumC, ModSizeJ, ModSizeC, RelocCostJ, RelocCostC
$load CapI
$load ManCost, RemanCost, ColCost, DispCost
$gdxIn