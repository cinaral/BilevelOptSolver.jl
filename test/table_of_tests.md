|BASBLib|BOLIB|notes|
|---|---|---|
|LP-LP|
|[mb_2007_01](./lp-lp/mb_2007_01.jl)||simplest|
|[mb_2007_02](./lp-lp/mb_2007_02.jl)||No solution|
|[bf_1982_02](./lp-lp/bf_1982_02.jl)|BardFalk1982Ex2|BASBLib version has upper bounds on x. our method fails without random initialization. doesn't work with smaller tolerances tol$<10^{-5}$|
|[ct_1982_01](./lp-lp/ct_1982_01.jl)|CandlerTownsley1982|BASBLib version has upper bounds on x and zero part of y is removed. our method requires forcing solution sets otherwise converges to false solution|
|LP-QP|
|[mb_2006_01](./lp-qp/mb_2006_01.jl)||non-convex f with 1 stationary point|
|[mb_2007_03](./lp-qp/mb_2007_03.jl)|||
|[mb_2007_04](./lp-qp/mb_2007_04.jl)|||
|[b_1991_02](./lp-qp/b_1991_02.jl)|Bard1991Ex1|BASBLib version has upper bounds on y|
|[as_1984_01](./lp-qp/as_1984_01.jl)|AiyoshiShimizu1984Ex2|our method requires forcing solution sets for some inits otherwise converges to false solution|
|QP-QP|
|[y_1996_02](./qp-qp/y_1996_02.jl)|Yezza1996Ex31|
|[d_1978_01](./qp-qp/d_1978_01.jl)|DeSilva1978|BASBLib version has upper bounds on x|
|[as_1981_01](./qp-qp/as_1981_01.jl)|||
|LP-NLP|||
|[mb_2007_05](./lp-nlp/mb_2007_05.jl)|MitsosBarton2006Ex35||
|[mb_2007_06](./lp-nlp/mb_2007_06.jl)|||
|[mb_2007_13](./lp-nlp/mb_2007_13.jl)|KleniatiAdjiman2014Ex3||
|[mb_2007_13v](./lp-nlp/mb_2007_13v.jl)|||
|[ka_2014_01](./lp-nlp/ka_2014_01.jl)||our method fails for x<0|
|NLP-NLP|
|[ka_2014_02](./nlp-nlp/ka_2014_02.jl)|KleniatiAdjiman2014Ex4||

|no optimal solutions|
|---|
|Zlobec2001b|
|MitsosBartonEx32|
|AnEtal|