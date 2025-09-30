using CSV

init_solver = "IPOPT"
solver = "IPOPT"

# you could simply use these instead, but this example is in the paper table order...
#include("../benchmarks/BASBLib_benchmark.jl")
#res = benchmark_BASBLib(; init_solver, solver);
#include("../benchmarks/BOLIB_benchmark.jl")
#res = benchmark_BOLIB(; init_solver, solver;

function split_list(list)
    items = split(strip(list), '\n')
end

function convert_to_id_list(merged_list, examples)
    ids = []
    list = split_list(merged_list)
    for name in list
        push!(ids, findfirst(isequal(String(name)), examples))
    end
    ids
end

include("../benchmarks/BOLIB_benchmark.jl")
LPLP_BOLIB_p1 = "MitsosBarton2006Ex31
MitsosBarton2006Ex32
ClarkWesterberg1988
LiuHart1994
Bard1984a
AnandalinghamWhite1990
Bard1991Ex2
ClarkWesterberg1990b
BardFalk1982Ex2
CandlerTownsley1982"
LPLP_BOLIB_p2 = "MershaDempe2006Ex2
WangJiaoLi2005
BenAyedBlair1990b
GlackinEtal2009
HaurieSavardWhite1990
Bard1984b
BenAyedBlair1990a
MershaDempe2006Ex1
TuyEtal1994
VisweswaranEtal1996
FrankeEtal2018Ex511
TuyEtal1993
BialasKarwan1984b
HuHuangZhang2009
LanWenShihLee2007
BialasKarwan1984a
TuyEtal2007Ex3"
LPQP_BOLIB_p1 = "MitsosBarton2006Ex34
TuyEtal2007
Bard1991Ex1
AiyoshiShimizu1984Ex2"
LPQP_BOLIB_p2 = "DempeEtal2012
MitsosBarton2006Ex33
FrankeEtal2018Ex521
LamparielloSagratella2017Ex23
HatzEtal2013
FrankeEtal2018Ex513
FrankeEtal2018Ex53
OutrataCervinka2009
DempeFranke2014Ex38
FloudasEtal2013"
QPQP_BOLIB_p1 = "Dempe1992b
LucchettiEtal1987
Yezza1996Ex31
ClarkWesterberg1990a
ShimizuAiyoshi1981Ex1
Bard1988Ex1
DempeDutta2012Ex31
DeSilva1978
FalkLiu1995
ShimizuAiyoshi1981Ex2
Bard1988Ex3
Bard1988Ex2"
QPQP_BOLIB_p2 = "HenrionSurowiec2011
LamparielloSagratella2017Ex32
MacalHurter1997
HendersonQuandt1958
GumusFloudas2001Ex4
SahinCiric1998Ex2
Yezza1996Ex41
CalamaiVicente1994a
AllendeStill2013
Dempe1992a
MuuQuy2003Ex1
Outrata1990Ex2a
Outrata1990Ex2b
Outrata1990Ex2c
Outrata1990Ex2d
Outrata1990Ex2e
Outrata1993Ex31
Outrata1993Ex32
Outrata1994Ex31
AnEtal2009
Outrata1990Ex1a
Outrata1990Ex1b
Outrata1990Ex1c
Outrata1990Ex1d
Outrata1990Ex1e
DempeFranke2011Ex41
DempeFranke2011Ex42
DempeLohse2011Ex31a
MuuQuy2003Ex2
BardBook1998
DempeLohse2011Ex31b
CalamaiVicente1994b
CalamaiVicente1994c
GumusFloudas2001Ex1
DesignCentringP4"
LPNLP_BOLIB_p1 = "MitsosBarton2006Ex35
MitsosBarton2006Ex36
MitsosBarton2006Ex310
MitsosBarton2006Ex311
MitsosBarton2006Ex313
MitsosBarton2006Ex315
MitsosBarton2006Ex316
KleniatiAdjiman2014Ex3
MitsosBarton2006Ex39
NieWangYe2017Ex34"
LPNLP_BOLIB_p2 = "MorganPatrone2006a
PaulaviciusAdjiman2017b
Zlobec2001b"
QPNLP_BOLIB_p1 = "MitsosBarton2006Ex322
DempeDutta2012Ex24
YeZhu2010Ex42
MitsosBarton2006Ex312
MitsosBarton2006Ex314
MitsosBarton2006Ex317
MitsosBarton2006Ex318
MitsosBarton2006Ex319
MitsosBarton2006Ex320
MitsosBarton2006Ex321
MitsosBarton2006Ex324
ShimizuEtal1997b
MitsosBarton2006Ex38
Colson2002BIPA2
Colson2002BIPA4"
QPNLP_BOLIB_p2 = "Mirrlees1999
YeZhu2010Ex43
Vogel2012
PaulaviciusAdjiman2017a
MitsosBarton2006Ex323
LuDebSinha2016e
DesignCentringP1
DesignCentringP2"
others_BOLIB = "Colson2002BIPA3
Colson2002BIPA1
NieWangYe2017Ex52
NieWangYe2017Ex57
Colson2002BIPA5
FloudasZlobec1998
NieWangYe2017Ex54
NieWangYe2017Ex58
MitsosBarton2006Ex326
KleniatiAdjiman2014Ex4
NieWangYe2017Ex61
LuDebSinha2016d
SinhaMaloDeb2014TP7
SinhaMaloDeb2014TP8
TollSettingP5
TollSettingP4
MitsosBarton2006Ex325
DesignCentringP3
NetworkDesignP1
NetworkDesignP2
MitsosBarton2006Ex327
MitsosBarton2006Ex328
TollSettingP1
RobustPortfolioP1
RobustPortfolioP2
TollSettingP2
TollSettingP3
LamparielloSagratella2017Ex31
LamparielloSagratella2017Ex35
LamparielloSagratella2017Ex33
Zlobec2001a
IshizukaAiyoshi1992a
CalveteGale1999P1
WanWangLv2011
LuDebSinha2016f"
example_ids =
    [
		convert_to_id_list(LPLP_BOLIB_p1, BOLIB.examples)
        convert_to_id_list(LPLP_BOLIB_p2, BOLIB.examples)
        convert_to_id_list(LPQP_BOLIB_p1, BOLIB.examples)
        convert_to_id_list(LPQP_BOLIB_p2, BOLIB.examples)
		convert_to_id_list(QPQP_BOLIB_p1, BOLIB.examples)
		convert_to_id_list(QPQP_BOLIB_p2, BOLIB.examples)
		convert_to_id_list(LPNLP_BOLIB_p1, BOLIB.examples)
		convert_to_id_list(LPNLP_BOLIB_p2, BOLIB.examples)
		convert_to_id_list(QPNLP_BOLIB_p1, BOLIB.examples)
		convert_to_id_list(QPNLP_BOLIB_p2, BOLIB.examples)
		convert_to_id_list(others_BOLIB, BOLIB.examples)
    ]
BOLIB_res = benchmark_BOLIB(; example_ids, init_solver, solver, do_force_dry_run=true);

include("../benchmarks/BASBLib_benchmark.jl")
LPLP_BASBLib = "as_2013_01
sib_1997_02
sib_1997_02v
b_1991_01v
s_1989_01
ct_1982_01"
LPQP_BASBLib = "mb_2007_03
mb_2006_01"
QPQP_BASBlib = "b_1998_04
b_1998_05
sc_1998_01
d_2000_01
b_1998_02
b_1998_03
b_1998_07
as_1981_01"
LPNLP_BASBLib = "mb_2007_13v
gf_2001_01
cg_1999_01"
QPNLP_BASBLib = "mb_2007_18v
mb_2007_22"
example_ids =
    [
        convert_to_id_list(LPLP_BASBLib, BASBLib.examples)
        convert_to_id_list(LPQP_BASBLib, BASBLib.examples)
        convert_to_id_list(QPQP_BASBlib, BASBLib.examples)
        convert_to_id_list(LPNLP_BASBLib, BASBLib.examples)
        convert_to_id_list(QPNLP_BASBLib, BASBLib.examples)
    ]
BASBLib_res = benchmark_BASBLib(; example_ids, init_solver, solver, do_force_dry_run=true);

res = [
    BOLIB_res
    BASBLib_res
]

CSV.write("BOLIB_BASBLib_ipopt.csv", res)
