# this code is building on the experimental code of Zdzislaw Jackiewicz in
# his matlab code dim18.m

function dim18(f, xspan, y0; abstol = 1e-6, reltol = 1e-6, pmax = 8)
    # Define fac, facmax and facmin used to compute new stepsize
    fac = 0.9
    facmax = 2.0
    facmin = 0.5
    # Define rmin and rmax used to make the decision about changing the order
    rmin = 0.8
    rmax = 1.2
    # Define the minimal number of rejections to reduce order by one
    nrmin = 3

    # Define coefficients of DIMSIMs in Nordsieck form
    ep = zeros(7, 7)
    # p=q=s=1
    c1 = Float64[0.0]
    A1 = Float64[0.0]
    P1 = Float64[1, 0]
    G1 = Float64[1; 1]
    Q1 = Float64[1 0; 0 0]
    va1 = [1/2, 1.0]
    fi1 = va1 - G1*(c1^1)/1
    EC[1] = Q1[1, :]*fi1
    ep[1, 1] = 1
    # p=q=s=2
    c2 = Float64[0, 1]
    A2 = Float64[0 0; 2 0]
    P2 = Float64[1 0 0; 1 -1 1/2]
    G2 = Float64[5/4, 1/4; 0, 1; -1, 1]
    Q2 = Float64[1 -1/2 1/4; 0 0 0; 0 0 0]
    va2 = Float64[1/6, 1/2, 1]
    fi2 = va2 - G2*(c2.^2)/2
    EC[2] = Q2[1, :]*fi2
    ep[2, 1] = 0.0
    ep[2, 2] = 1/2
    # p=q=s=3
    c3 = Float64[0, 1/2, 1]
    A3 = Float64[0 0 0; 1 0 0; 1/4, 1, 0]
    P3 = Float64[1 0 0 0; 1 -1/2 1/8 1/48; 1 -1/4 0 1/24]
    G3 = Float64[5/4 1/3 1/6; 0 0 1; 1 -4 3; 4 -8 4]
    Q3 = Float64[1 -3/4 1/6 1/24; 0 0 0 0; 0 0 0 0; 0 0 0 0]
    va3 = Float64[1/24, 1/6, 1/2, 1]
    fi3 = va3 - G3*(c3.^3)/6
    EC[3] = Q3[1, :]*fi3
    ep[3, 1] = 0
    ep[3, 2] = 1/12
    ep[3, 3] = 1/2
    # p=q=s=4
    c4 = Float64[0, 1/3, 2/3, 1]
    A4 = Float64[
         0 0 0 0;
         0.3739348246366938456769753 0 0 0;
         0.2949848976674274114681988 0.4816828233453680326779301 0 0;
        -0.6903740089436590462138745 2.271260297677338864939909 -0.2257249931856950422183244 0
    ]
    P4 = Float64[
        1 0 0 0 0;
        1 -0.04060149130336049 0.05555555555555555 0.006172839506172839 0.00051440329218107;
        1 -0.1100010543461287 0.06166128110709955 0.02262255919686226 0.005257101913505959;
        1 -0.3551612955479846 -0.1066034371019828 0.0906466486147468 0.03879345461610014
    ]
    G4 = Float64[
        2.944452437155611 -6.16786637939422 4.995564740173299 -1.024830488098152;
        0 0 0 1;
        -1 4.5 -9 5.5;
        -9 36 -45 18;
        -27 81 -81 27
    ]
    Q4=[1,0.2526796901634609, 0.2504094544473587, ...
        -0.0883843438008673,0.003850442201592136; ...
        0,0,0,0,0; ...
        0,0,0,0,0; ...
        0,0,0,0,0; ...
        0,0,0,0,0];
    va4=[1/120;1/24;1/6;1/2;1];
    fi4=va4-G4*(c4.^4)/24;
    EC(4)=Q4(1,:)*fi4;
    ep(4,1)=0;
    ep(4,2)=1/108;
    ep(4,3)=11/108;
    ep(4,4)=1/2;

    #    p=q=s=5
    c5=[0;1/4;2/4;3/4;1];
    A5=[0,0,0,0,0; ...
        1.176528170310668,0,0,0,0; ...
        1.980579346323319, 0.4017181027378084, 0, 0, 0; ...
        3.053283510839239, 0.7349961462028882, 0.2672357626475791, 0, 0; ...
        1.419332526946769, 2.653489747312533, -2.277853294546826, ...
        1.197890508817277, 0];
    P5=[1,0,0,0,0,0; ...
        1, -0.926528170310669, 0.03125, ...
        0.002604166666666666, 0.0001627604166666666, 8.13802083333333*10^(-6); ...
        1, -1.882297449061127, 0.02457047431554787, 0.00827964262277682, ...
        0.00155802577412029, 0.0001950328608825181; ...
        1, -3.305515419689706, -0.03611691787451161, 0.01393940010021235, ...
        0.005702129564105414, 0.001161984318267552; ...
        1, -1.992859488529754, 0.07713632883232156, 0.0315700682766439, ...
        -0.002014862315115688, -0.001959141967572073];
    G5=[3.163023914364555, 1.97433924605054, -0.810425120055628, ...
        0.5409220188020185, 0.05507865419631682; ...
        0,0,0,0,1; ...
        1,-16/3,12,-16,25/3; ...
        44/3,-224/3,152,-416/3,140/3; ...
        96, -448, 768, -576, 160; ...
        256, -1024, 1536, -1024, 256];
    Q5=[1, -3.922938713357803, -0.0491424197826518, 0.02659806034831467, ...
        0.006195659734832572, 0.0006962310673207849; ...
        0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0];
    va5=[1/720;1/120;1/24;1/6;1/2;1];
    fi5=va5-G5*(c5.^5)/120;
    EC(5)=Q5(1,:)*fi5;
    ep(5,1)=0;
    ep(5,2)=1/1280;
    ep(5,3)=5/384;
    ep(5,4)=7/64;
    ep(5,5)=1/2;
    #    p=q=s=6
    c6=[0;1/5;2/5;3/5;4/5;1];
    A6=[0,0,0,0,0,0; ...
        -0.3787322586975567,0,0,0,0,0; ...
        3.415474814317656, 1.075852942341251, 0, 0, 0, 0; ...
        3.517453259414935, 0.493910093894091, 0.4832074216261553, 0, 0, 0; ...
        -0.3980974661620062, -1.988292762808572, 1.043842001373579, ...
        0.253273266522039, 0, 0; ...
        1.348865151959959, -2.741610520631099, 1.388314181048896, ...
        0.03275797313010722, 0.3645870239782025, 0];
    P6=[1,0,0,0,0,0,0; ...
        1,287342184/496502795,1/50,1/750,1/15000,1/375000,1/11250000; ...
        1, -4.091327756658908, -0.1351705884682503, -0.01085039218015837, ...
        -0.0003678039231216692, 0.00001360980384391653, 2.819947709312216*10^(-6); ...
        1, -3.894570774935182, -0.1120649874292803, -0.01253479560797424, ...
        -0.000412759289204445, 0.0000996514106724949, 0.00002224920643751716; ...
        1, 1.88927496107496, 0.1481577920990593, -0.003997359494348594, ...
        -0.0005344285923668219, 0.0003821124101697426, 0.0001161954087662182; ...
        1, 0.6070861905139336, 0.1816720286460348, -0.01213054024106717, ...
        -0.001777250315837238, 0.0006360606441260151, 0.0002609375734296115];
    G6=[5.260629843769215, -1.616594220242467, 1.675267728346385, ...
        -0.6076770669598051, 0.736220098847675, -0.07638264443638971; ...
        0,0,0,0,0,1; ...
        -1,25/4,-50/3,25,-25,137/12; ...
        -125/6,1525/12,-325,2675/6,-1925/6,375/4; ...
        -875/4,5125/4,-6125/2,7375/2,-8875/4,2125/4; ...
        -1250,6875,-15000,16250,-8750,1875; ...
        -3125,15625,-31250,31250,-15625,3125];
    Q6=[1, -4.371463739324614, 0.005224558244071898, ...
        -0.02304010455649106, -0.002264696760421819, 0.0005533973638788859, ...
        0.0002701687204304608; ...
        0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0];
    va6=[1/5040;1/720;1/120;1/24;1/6;1/2;1];
    fi6=va6-G6*(c6.^6)/720;
    EC(6)=Q6(1,:)*fi6;
    ep(6,1)=0;
    ep(6,2)=1/18750;
    ep(6,3)=137/112500;
    ep(6,4)=3/200;
    ep(6,5)=17/150;
    ep(6,6)=1/2;

    #    p=q=s=7
    c7=[0;1/6;2/6;3/6;4/6;5/6;1];
    A7=[0,0,0,0,0,0,0; ...
        -0.12432301443470899101, 0, 0, 0, 0, 0, 0; ...
        -0.20080968439550891353, 0.4391569966420439819, 0, 0, 0, 0, 0; ...
        0.08616601501555395074, -0.04036381021150520290, 0.4215108431387438460, ...
        0, 0, 0, 0; ...
        0.5012952090985448529, -0.4785488252539539821, ...
        0.3811857626426213823, 0.3273204640053006166, 0, 0, 0; ...
        0.7964645243641904650, -0.8959743476016479602, 0.7564909788229435454, ...
        -0.09019470959572305346, 0.3872789952897844843, 0, 0; ...
        0.5942581855169782141, 0.15709815967798328794, 0.08129751241083920929, ...
        0.11216039145279535853, 0.15567633255920378469, 0.3384770251065685742, 0];
    P7=[1,0,0,0,0,0,0,0; ...
        1, 0.29098968110137565768, 0.013888888888888888889, ...
        0.0007716049382716049383, 0.00003215020576131687243, ...
        1.0716735253772290809*10^(-6), 2.9768709038256363359*10^(-8), ...
        7.087787866251515086*10^(-10); ...
        1, 0.09498602108679826493, ...
        -0.017637277218118441433, 0.00007343677503333975714, ...
        0.00017554758489554219501, 0.000020174565008507673763, ...
        1.4345644516629520274*10^(-6), 7.765054783286785900*10^(-8); ...
        1, 0.03268695205720740615, ...
        -0.008776312677663748187, -0.0020233272547704191825,  ...
        0.00003339279714626113798, 0.00004488780606968890389, ...
        7.289541355636683525*10^(-6), 7.482394315355554629*10^(-7); ...
        1, -0.06458594382584620305, ...
        0.011258206881023783490, -0.006062817302898354738, ...
        -0.0005724715566076134710, 0.00006429888026864696670, ...
        0.000024133562076881641897, 3.797334624446901293*10^(-6); ...
        1, -0.12073210794621414771, ...
        0.031298978486187570983, -0.007920231410010083953, ...
        -0.0011303127745192270740, 0.00003604468994887979844, ...
        0.00003864418730502681580, 8.694153244262549844*10^(-6); ...
        1, -0.4389676067243684288, 0.004789197562049034676, ...
        -0.006173315013999628643, -0.0016271005017749117589, ...
        -0.00008822440498615491390, 0.00005233318150022978818, ...
        0.000019399196523704375274];
    G7=[-0.25833088744103508910, 2.6466968223136745127, ...
        -2.1476803074095682625, 1.9316924755337628148, -1.4976252916292613939, ...
        1.2618657088457168638, -0.20093642140062518816; ...
        0,0,0,0,0,0,1; ...
        1, -36/5, 45/2, -40, 45, -36, 147/10; ...
        137/5, -972/5, 594, -1016, 1053, -3132/5, 812/5; ...
        405, -2808, 8289, -13392, 12447, -6264, 1323; ...
        3672, -24624, 69336, -104544, 88776, -40176, 7560; ...
        19440, -124416, 330480, -466560, 369360, -155520, 27216; ...
        46656, -279936, 699840, -933120, 699840, -279936, 46656];
    Q7=[1, -0.7356820988126642575, ...
        -0.04327041390060234104, 0.0028870172620625095294, ...
        -0.0016233635195644772625, -0.0003347252333547855415, ...
        0.00004865116250519555581, 0.00003525374418906881682; ...
        0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0];
    va7=[1/40320;1/5040;1/720;1/120;1/24;1/6;1/2;1];
    fi7=va7-G7*(c7.^7)/5040
    EC[7] = Q7[1, :]*fi7
    ep[7, 1] = 0;
    ep[7, 2] = 1/326592
    ep[7, 3] = 7/77760
    ep[7, 4] = 29/19440
    ep[7, 5] = 7/432
    ep[7, 6] = 25/216
    ep[7, 7] = 1/2

    #    p=q=s=8
    c8=[0;1/7;2/7;3/7;4/7;5/7;6/7;1];
    A8=[0,0,0,0,0,0,0,0; ...
        1.2860282000085701934, 0, 0, 0, 0, 0, 0, 0; ...
        -8.107976308060821188, 0.5286183572282346903, 0, 0, 0, 0, 0, 0; ...
        -13.936070304367901018, -0.12031917097591100229, 0.3947575778043359273, ...
        0, 0, 0, 0, 0; ...
        -14.666442638286900761, -1.3340533277983598392, ...
        0.4014562781039848100, 0.29602776935058129555, 0, 0, 0, 0; ...
        -15.366227398890298965, -3.418800902216350148, 0.6550270118416955847, ...
        0.19414843372705584818, 0.26968528159293970545, 0, 0, 0; ...
        -21.554460383192100535, -5.703764562751610121, 0.5295796362522783580, ...
        0.7534525920384679702, -0.14815838628387403519, 0.3255280507569289984, ...
        0, 0; ...
        -29.080046065331000938, -4.964275870797740592, ...
        -2.2828012028605800906, 3.533190437782849788, -1.3936788059280600252,  ...
        0.3489546294105012917, 0.29594458482908110432, 0];
    P8=[1,0,0,0,0,0,0,0,0; ...
        1, -1.1431710571514273363, 0.010204081632653061224, ...
        0.0004859086491739552964, 0.000017353880327641260586, ...
        4.958251522183217310*10^(-7), 1.1805360767102898358*10^(-8), ...
        2.4092572994087547669*10^(-10), 4.302245177515633512*10^(-12); ...
        1, 7.865072236546872212, -0.03470058164484985372, ...
        -0.0015067956762842218160, 0.000020801853352933303316, ...
        6.692825160653193034*10^(-6), 4.934408116564968561*10^(-7), ...
        2.4597963017239474377*10^(-8), 9.740170018686432681*10^(-10); ...
        1, 14.090203326110904664, -0.003762691682231142154, ...
        -0.0017652740257872142087, -0.00007040053867912223451, ...
        0.000012964304264953474652, 2.4023817151484592099*10^(-6), ...
        2.3006862282347988808*10^(-7), 1.6082289623251521997*10^(-8); ...
        1, 15.874440490059265923, 0.11227351519939845568, ...
        0.0011387483555738424322, -0.0003534934531788487675, ...
        3.291101240461874957*10^(-6), 6.979489696797522317*10^(-6), ...
        1.1121116698747153461*10^(-6), 1.1391465600327489778*10^(-7); ...
        1, 18.380453288230672260, 0.3190392480991883324, ...
        0.007028299981632472138, -0.0009727173730161769440, ...
        -0.00004410266939505382812, 0.000015442944556953841532, ...
        3.656352861934463038*10^(-6), 4.943544196008848963*10^(-7); ...
        1, 26.654965910322766507, 0.5600956252730102960,  ...
        0.013494080476829900765, -0.0018460873902036699560, ...
        -0.00012415357898412204540, 0.000025269435420637103162, ...
        7.744367157404060043*10^(-6), 1.2718005850811728574*10^(-6); ...
        1, 34.54271229289494946, 0.6406545044079426323, ...
        0.015827450829257699507, -0.0023165071654290093980, ...
        -0.00016240422994434372568, 0.000027763279873144490238, ...
        9.808713943675449713*10^(-6), 1.9850738649671162216*10^(-6)];
    G8=[-29.346017024479626510, -5.262452999422931031, ...
        -2.1025831257271133070, 3.552139013365491294, -1.5839739310677373933, ...
        0.5315770849396886481, 0.22029175309245071925, 0.011983723374869140579; ...
        0,0,0,0,0,0,0,1; ...
        -1, 49/6, -147/5, 245/4, -245/3, 147/2, -49, 363/20; ...
        -343/10, 49931/180, -9849/10, 2009, -46501/18, 43071/20, -10927/10, 22981/90; ...
        -9947/15, 634207/120, -91924/5, 872935/24, -133427/3, 1347647/40, ...
        -218834/15, 331681/120; ...
        -16807/2, 196882/3, -444185/2, 422576, -2926819/6, 340942, ...
        -266511/2, 67228/3; ...
        -420175/6, 1596665/3, -3479049/2, 9495955/3, ...
        -20756645/6, 2268945, -4958065/6, 386561/3; ...
        -352947, 2588278, -8117781, 14117880, -14706125, 9176622, -3176523, ...
        470596; ...
        -823543, 5764801, -17294403, 28824005, -28824005, 17294403, ...
        -5764801, 823543];
    Q8=[1, 34.97903550592490844, 0.6547965611677716853, ...
        0.016054264615000655678, -0.0023522231148220533850, ...
        -0.00016709809517544116894, 0.000028239128344368648942, ...
        0.000010053192572845828464, 2.0079727295695647434*10^(-6); ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0, 0, 0, 0, 0];
    va8=[1/362880;1/40320;1/5040;1/720;1/120;1/24;1/6;1/2;1];
    fi8 = va8 - G8*(c8.^8)/40320
    EC[8] = Q8[1, :]*fi8

    # Define the interval of integration [a,b]
    a, b = xspan
    # Define the dimension of the system of equations m
    m = length(y0)
    # Define the absolute and relative error tolerances
    for j = 1:m
        Atol[j] = abstol;
        Rtol[j] = reltol;
    end
    # Define the initial gridpoint x0=x(0)
    x0 = a

    ystart = y0

    nfe = 0

    # Initial stepsize using the approach described in the paper by
    # Butcher, Chartier and Jackiewicz, Experiments with a variable-
    # order type 1 DIMSIM code, Numerical Algorithms 22(1999), 237-261.
    # 1. Define sc
    for j = 1:m
        sc[j] = Atol[j] + abs(y0[j])*Rtol[j]
    end
    # 2. Compute d0 and d1
    sum = 0.0
    for j = 1:m
        sum += (y0[j]/sc[j])^2
    end
    d0 = sqrt(sum/m)
    f0 = f(x0, y0)
    sum = 0
    for j = 1:m
        sum += (f0[j]/sc[j])^2
    end
    d1 = sqrt(sum/m)
    # 3. Compute d2
    h0 = 1/d1
    d2 = (f(x0 + h0, y0 + h0*f0) - f0)/h0
    nfe += 2
    sum = 0.0
    for j = 1:m
        sum += (d2[j]/sc[j])^2
    end
    sum = sqrt(sum/m)
    h0 = sqrt(2/sum)
    h0 = min(0.1*(b - a), h0)


    # Compute starting vector z0
    z0 = [y0, h0*f0]

    # Initialize the algorithm
    xold = x0
    hol = h0
    hnew = h0
    pold = 1
    pnew = 1
    xnew = xold + hnew
    del = hnew/hol
    Yold = y0 - hol*f0
    Fold = f(xold - hol, Yold)
    zold = z0;
    # Define modified starting vector
    zold = [y0, h0*Fold]

    ns = 0             # Number of steps
    nr = 0             # Number of rejected steps
    nrl = 0            # Number of local rejections
    xout = a           # Initial output of independent variable
    yout = ystart      # Initial condition
    h = []             # Stepsizes
    p = []             # Orders
    xr = []            # Values of independent variables where rejected steps occur
    hr = []            # Rejected stepsizes
    pr = []            # Rejected orders
    esth = []          # Estimates of local discretization errors
    nesth = []         # Norm of esth

    ############################################################################
    # Main algorithm
    ############################################################################
    Im = sparse(1:m, 1:m, ones(m, 1), m, m, m)

    while xold < b
        # Step of order 1
        if pnew == 1
            zold = zold(1:(pnew+1)*m);
            D = diag(del.^(0:pnew));
            Ynew = kron(P1*D,Im)*zold;
            Fnew = feval(f,xold+c1(1)*hnew,Ynew,varargin{:});
            nfe = nfe + 1;
            znew = hnew*kron(G1,Im)*Fnew+kron(Q1*D,Im)*zold;
            # Compute local error estimate
            bet1z = [del/2];
            gam1z = -del^2/2;
            est = hnew*kron(bet1z', Im)*Fnew + gam1z*zold((pnew*m + 1):(pnew+1)*m)
        end

        # Step of order 2
        if pnew == 2
            # Define vector zold
            if (pnew-pold==1)
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if (pnew-pold==-1)
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P2*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c2(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A2(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c2(2)*hnew,Ynew2,varargin{:});
            nfe=nfe+2;
            Ynew=[Ynew1;Ynew2];
            Fnew=[Fnew1;Fnew2];
            znew=hnew*kron(G2,Im)*Fnew+kron(Q2*D,Im)*zold;
            # Compute local error estimate
            bet2z=[-(3+del)/(12*(1+del));(3+del)/(12*(1+del))];
            gam2z=-del^2*(3+del)/(12*(1+del));
            est=hnew*kron(bet2z',Im)*Fnew+gam2z*zold(pnew*m+1:(pnew+1)*m);;
        end

        # Step of order 3

        if pnew == 3
            # Define vector zold
            if pnew - pold == 1
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if pnew - pold == -1
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P3*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c3(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A3(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c3(2)*hnew,Ynew2,varargin{:});
            Ynew3=hnew*(A3(3,1)*Fnew1+A3(3,2)*Fnew2)+ztemp(2*m+1:3*m);
            Fnew3=feval(f,xold+c3(3)*hnew,Ynew3,varargin{:});
            nfe=nfe+3;
            Ynew=[Ynew1;Ynew2;Ynew3];
            Fnew=[Fnew1;Fnew2;Fnew3];
            znew=hnew*kron(G3,Im)*Fnew+kron(Q3*D,Im)*zold;
            # Compute local error estimate
            bet3z=[(2+del)/(18*del);-(2+del)/(9*del);(2+del)/(18*del)];
            gam3z=-del^2*(2+del)/72;
            est=hnew*kron(bet3z',Im)*Fnew+gam3z*zold(pnew*m+1:(pnew+1)*m);
        end

        # Step of order 4
        if (pnew==4)
            # Define vector zold
            if (pnew-pold==1)
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if (pnew-pold==-1)
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P4*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c4(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A4(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c4(2)*hnew,Ynew2,varargin{:});
            Ynew3=hnew*(A4(3,1)*Fnew1+A4(3,2)*Fnew2)+ztemp(2*m+1:3*m);
            Fnew3=feval(f,xold+c4(3)*hnew,Ynew3,varargin{:});
            Ynew4=hnew*(A4(4,1)*Fnew1+A4(4,2)*Fnew2+A4(4,3)*Fnew3)+ztemp(3*m+1:4*m);
            Fnew4=feval(f,xold+c4(4)*hnew,Ynew4,varargin{:});
            nfe=nfe+4;
            Ynew=[Ynew1;Ynew2;Ynew3;Ynew4];
            Fnew=[Fnew1;Fnew2;Fnew3;Fnew4];
            znew=hnew*kron(G4,Im)*Fnew+kron(Q4*D,Im)*zold;
            # Compute local error estimate
            bet4z=[-0.45*del/(1+del);1.35*del/(1+del);-1.35*del/(1+del);0.45*del/(1+del)];
            gam4z=-0.01666666666666669*del^5/(1+del);
            est=hnew*kron(bet4z',Im)*Fnew+gam4z*zold(pnew*m+1:(pnew+1)*m);
        end

        # Step of order 5
        if pnew == 5
            # Define vector zold
            if pnew - pold == 1
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if pnew - pold == -1
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P5*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c5(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A5(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c5(2)*hnew,Ynew2,varargin{:});
            Ynew3=hnew*(A5(3,1)*Fnew1+A5(3,2)*Fnew2)+ztemp(2*m+1:3*m);
            Fnew3=feval(f,xold+c5(3)*hnew,Ynew3,varargin{:});
            Ynew4=hnew*(A5(4,1)*Fnew1+A5(4,2)*Fnew2+A5(4,3)*Fnew3)+ztemp(3*m+1:4*m);
            Fnew4=feval(f,xold+c5(4)*hnew,Ynew4,varargin{:});
            Ynew5=hnew*(A5(5,1)*Fnew1+A5(5,2)*Fnew2+A5(5,3)*Fnew3+A5(5,4)*Fnew4)+ ...
                  ztemp(4*m+1:5*m);
            Fnew5=feval(f,xold+c5(5)*hnew,Ynew5,varargin{:});
            nfe=nfe+5;
            Ynew=[Ynew1;Ynew2;Ynew3;Ynew4;Ynew5];
            Fnew=[Fnew1;Fnew2;Fnew3;Fnew4;Fnew5];
            znew=hnew*kron(G5,Im)*Fnew+kron(Q5*D,Im)*zold;
            # Compute local error estimate
            bet5z=[0.7111111111124374*del/(1+del); ...
            -2.84444444444975*del/(1+del); ...
            4.266666666674624*del/(1+del); ...
            -2.84444444444975*del/(1+del); ...
            0.7111111111124374*del/(1+del)];
            gam5z=-0.002777777777782958*del^6/(1+del);
            est=hnew*kron(bet5z',Im)*Fnew+gam5z*zold(pnew*m+1:(pnew+1)*m);
        end

        # Step of order 6
        if pnew == 6
            # Define vector zold
            if pnew - pold == 1
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if pnew - pold == -1
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P6*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c6(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A6(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c6(2)*hnew,Ynew2,varargin{:});
            Ynew3=hnew*(A6(3,1)*Fnew1+A6(3,2)*Fnew2)+ztemp(2*m+1:3*m);
            Fnew3=feval(f,xold+c6(3)*hnew,Ynew3,varargin{:});
            Ynew4=hnew*(A6(4,1)*Fnew1+A6(4,2)*Fnew2+A6(4,3)*Fnew3)+ztemp(3*m+1:4*m);
            Fnew4=feval(f,xold+c6(4)*hnew,Ynew4,varargin{:});
            Ynew5=hnew*(A6(5,1)*Fnew1+A6(5,2)*Fnew2+A6(5,3)*Fnew3+A6(5,4)*Fnew4)+ ...
                  ztemp(4*m+1:5*m);
            Fnew5=feval(f,xold+c6(5)*hnew,Ynew5,varargin{:});
            Ynew6=hnew*(A6(6,1)*Fnew1+A6(6,2)*Fnew2+A6(6,3)*Fnew3+A6(6,4)*Fnew4 + ...
                  A6(6,5)*Fnew5)+ztemp(5*m+1:6*m);
            Fnew6=feval(f,xold+c6(6)*hnew,Ynew6,varargin{:});
            nfe=nfe+6;
            Ynew=[Ynew1;Ynew2;Ynew3;Ynew4;Ynew5;Ynew6];
            Fnew=[Fnew1;Fnew2;Fnew3;Fnew4;Fnew5;Fnew6];
            znew=hnew*kron(G6,Im)*Fnew+kron(Q6*D,Im)*zold;
            # Compute local error estimate
            bet6z=[-1.265588115871886*del/(1+del); ...
            6.327940579359435*del/(1+del); ...
            -12.65588115871887*del/(1+del); ...
            12.65588115871887*del/(1+del); ...
            -6.327940579359435*del/(1+del); ...
            1.265588115871886*del/(1+del)];
            gam6z=-0.0004049881970790037*del^7/(1+del);
            est=hnew*kron(bet6z',Im)*Fnew+gam6z*zold(pnew*m+1:(pnew+1)*m);
        end

        # Step of order 7
        if pnew == 7
            # Define vector zold
            if pnew - pold == 1
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if pnew - pold == -1
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P7*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c7(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A7(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c7(2)*hnew,Ynew2,varargin{:});
            Ynew3=hnew*(A7(3,1)*Fnew1+A7(3,2)*Fnew2)+ztemp(2*m+1:3*m);
            Fnew3=feval(f,xold+c7(3)*hnew,Ynew3,varargin{:});
            Ynew4=hnew*(A7(4,1)*Fnew1+A7(4,2)*Fnew2+A7(4,3)*Fnew3)+ztemp(3*m+1:4*m);
            Fnew4=feval(f,xold+c7(4)*hnew,Ynew4,varargin{:});
            Ynew5=hnew*(A7(5,1)*Fnew1+A7(5,2)*Fnew2+A7(5,3)*Fnew3+A7(5,4)*Fnew4)+ ...
                  ztemp(4*m+1:5*m);
            Fnew5=feval(f,xold+c7(5)*hnew,Ynew5,varargin{:});
            Ynew6=hnew*(A7(6,1)*Fnew1+A7(6,2)*Fnew2+A7(6,3)*Fnew3+A7(6,4)*Fnew4 + ...
                  A7(6,5)*Fnew5)+ztemp(5*m+1:6*m);
            Fnew6=feval(f,xold+c7(6)*hnew,Ynew6,varargin{:});
            Ynew7=hnew*(A7(7,1)*Fnew1+A7(7,2)*Fnew2+A7(7,3)*Fnew3+A7(7,4)*Fnew4 + ...
                  A7(7,5)*Fnew5+A7(7,6)*Fnew6)+ztemp(6*m+1:7*m);
            Fnew7=feval(f,xold+c7(7)*hnew,Ynew7,varargin{:});
            nfe=nfe+7;
            Ynew=[Ynew1;Ynew2;Ynew3;Ynew4;Ynew5;Ynew6;Ynew7];
            Fnew=[Fnew1;Fnew2;Fnew3;Fnew4;Fnew5;Fnew6;Fnew7];
            znew=hnew*kron(G7,Im)*Fnew+kron(Q7*D,Im)*zold;
            # Compute local error estimate
            bet7z=[2.325890162108168*del/(1+del); ...
            -13.95534097264901*del/(1+del); ...
            34.88835243162252*del/(1+del); ...
            -46.51780324216336*del/(1+del); ...
            34.88835243162252*del/(1+del); ...
            -13.95534097264901*del/(1+del); ...
            2.325890162108168*del/(1+del)];
            gam7z=-0.00004985189819333351*del^8/(1+del);
            est=hnew*kron(bet7z',Im)*Fnew+gam7z*zold(pnew*m+1:(pnew+1)*m);
        end

        # Step of order 8
        if pnew == 8
            # Define vector zold
            if pnew - pold == 1
                zold=[zold(1:pnew*m);est/EC(pnew-1)];
            end
            if pnew - pold == -1
                zold=zold(1:(pnew+1)*m);
            end
            D=diag(del.^(0:pnew));
            ztemp=kron(P8*D,Im)*zold;
            Ynew1=ztemp(1:m);
            Fnew1=feval(f,xold+c8(1)*hnew,Ynew1,varargin{:});
            Ynew2=hnew*A8(2,1)*Fnew1+ztemp(m+1:2*m);
            Fnew2=feval(f,xold+c8(2)*hnew,Ynew2,varargin{:});
            Ynew3=hnew*(A8(3,1)*Fnew1+A8(3,2)*Fnew2)+ztemp(2*m+1:3*m);
            Fnew3=feval(f,xold+c8(3)*hnew,Ynew3,varargin{:});
            Ynew4=hnew*(A8(4,1)*Fnew1+A8(4,2)*Fnew2+A8(4,3)*Fnew3)+ztemp(3*m+1:4*m);
            Fnew4=feval(f,xold+c8(4)*hnew,Ynew4,varargin{:});
            Ynew5=hnew*(A8(5,1)*Fnew1+A8(5,2)*Fnew2+A8(5,3)*Fnew3+A8(5,4)*Fnew4)+ ...
                  ztemp(4*m+1:5*m);
            Fnew5=feval(f,xold+c8(5)*hnew,Ynew5,varargin{:});
            Ynew6=hnew*(A8(6,1)*Fnew1+A8(6,2)*Fnew2+A8(6,3)*Fnew3+A8(6,4)*Fnew4 + ...
                  A8(6,5)*Fnew5)+ztemp(5*m+1:6*m);
            Fnew6=feval(f,xold+c8(6)*hnew,Ynew6,varargin{:});
            Ynew7=hnew*(A8(7,1)*Fnew1+A8(7,2)*Fnew2+A8(7,3)*Fnew3+A8(7,4)*Fnew4 + ...
                  A8(7,5)*Fnew5+A8(7,6)*Fnew6)+ztemp(6*m+1:7*m);
            Fnew7=feval(f,xold+c8(7)*hnew,Ynew7,varargin{:});
            Ynew8=hnew*(A8(8,1)*Fnew1+A8(8,2)*Fnew2+A8(8,3)*Fnew3+A8(8,4)*Fnew4 + ...
                  A8(8,5)*Fnew5+A8(8,6)*Fnew6+A8(8,7)*Fnew7)+ztemp(7*m+1:8*m);
            Fnew8=feval(f,xold+c8(8)*hnew,Ynew8,varargin{:});
            nfe=nfe+8;
            Ynew=[Ynew1;Ynew2;Ynew3;Ynew4;Ynew5;Ynew6;Ynew7;Ynew8];
            Fnew=[Fnew1;Fnew2;Fnew3;Fnew4;Fnew5;Fnew6;Fnew7;Fnew8];
            znew=hnew*kron(G8,Im)*Fnew+kron(Q8*D,Im)*zold;
            # Compute local error estimate
            bet8z=[-4.315167205300165*del/(1+del); ...
            30.20617043710115*del/(1+del); ...
            -90.6185113113035*del/(1+del); ...
            151.0308521855057*del/(1+del); ...
            -151.0308521855057*del/(1+del); ...
            90.6185113113035*del/(1+del); ...
            -30.20617043710115*del/(1+del); ...
            4.315167205300165*del/(1+del)];
            gam8z=-5.239759436119504*10^(-6)*del^9/(1+del);
            est=hnew*kron(bet8z',Im)*Fnew+gam8z*zold(pnew*m+1:(pnew+1)*m);
        end

        y0=zold(1:m);
        y1=znew(1:m);
        for j=1:m
            sc(j)=Atol(j)+max(abs(y0(j)),abs(y1(j)))*Rtol(j);
        end
        sum=0;
        for j=1:m
            sum=sum+(est(j)/sc(j))^2;
        end
        err=sqrt(sum/m);

        if err <= 1 # The step is accepted
            ns = ns + 1
            xout=[xout,xnew];
            yout=[yout,znew(1:m)];
            h=[h,hnew];
            p=[p,pnew];
            xold=xnew;
            if xnew < b
                hol=hnew;
                hnew=hnew*min(facmax,max(facmin,fac*(1/err)^(1/(pnew+1))));
                pold=pnew;
                if (pnew==1)
                    ratio=norm(est)/norm(znew(pnew*m+1:(pnew+1)*m));
                end
                if (pnew>1)
                    ratio=norm(est)/norm(znew(pnew*m+1:(pnew+1)*m)*EC(pnew-1));
                end
                if ratio < rmin && nrl == 0 && pnew < pmax
                    # Correct z^[n]
                    for j = 2:pnew+1
                        znew((j-1)*m+1:j*m)=znew((j-1)*m+1:j*m)+ep(pnew,j-1)*est/EC(pnew);
                    end
                    pnew=pnew+1;
                end
                if ratio > rmax && pnew > 1
                    pnew=pnew-1;
                end
                del=hnew/hol;
                xnew=xnew+hnew;
                zold=znew;
            end
            nrl=0;
        end

        if err > 1 # The step is rejected
            nrl=nrl+1;
            nr=nr+1;
            xr=[xr,xold];
            hr=[hr,hnew];
            pr=[pr,pnew];
            hnew=hnew*min(facmax,max(facmin,fac*(1/err)^(1/(pnew+1))));
            pold=pnew;
            if (nrl==nrmin && pnew>1)
                pnew=pnew-1;
            end
            if (nrl>nrmin)
                pnew=1;
            end
            del=hnew/hol;
            xnew=xold+hnew;
        end
    end

    return xout, yout
end
