#
#   Strang: r=2, s=1
#   SUZ90:  r=4, s=5  Suzuki
#   SS05: r=6, s=13 Sofroniu & Spaletta
#   SS05: r=8, s=21 Sofroniu & Spaletta
#   SS05: r=10, s=35 Sofroniu & Spaletta
#


function Splitting_Coefficients(r,h)

    Typeh=typeof(h)
    s=1

    if r==2
       # 2nd-order Strang 
       s=1
       a1=parse(BigFloat,"1.0")
       gamma=[a1]

    elseif r==4
        # 4th-order: Suzuki/Creutz&Gocksch (SUZ90)
        s=5
        a1=1/(4-4^(1/3))
        a2=a1
        a3=1-2*(a1+a2)
        a4=a2
        a5=a1
        gamma=[a1,a2,a3,a4,a5]

    elseif r==6
        # 6th-order: (SS05)
        s=13
        a1=parse(BigFloat,"0.13861930854051695245808013042625")
        a2=parse(BigFloat,"0.13346562851074760407046858832209")
        a3=parse(BigFloat,"0.13070531011449225190542755785015")
        a4=parse(BigFloat,"0.12961893756907034772505366537091")
        a5=parse(BigFloat,"-0.35000324893920896516170830911323")
        a6=parse(BigFloat,"0.11805530653002387170273438954049")
        a7=1-2*(a1+a2+a3+a4+a5+a6)  # =0.39907751534871587459988795520665
        a8=a6;    a9=a5;    a10=a4;    a11=a3;    a12=a2;    a13=a1
        gamma=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13]

    elseif r==8
        # 8th-order: (SS05)
        s=21
        a1=parse(BigFloat,"0.10647728984550031823931967854896")
        a2=parse(BigFloat,"0.10837408645835726397433410591546")
        a3=parse(BigFloat,"0.35337821052654342419534541324080")
        a4=parse(BigFloat,"-0.23341414023165082198780281128319")
        a5=parse(BigFloat,"-0.24445266791528841269462171413216")
        a6=parse(BigFloat,"0.11317848435755633314700952515599")
        a7=parse(BigFloat,"0.11892905625000350062692972283951")
        a8=parse(BigFloat,"0.12603912321825988140305670268365")
        a9=parse(BigFloat,"0.12581718736176041804392391641587")
        a10=parse(BigFloat,"0.11699135019217642180722881433533")
        a11=1-2*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10) # =-0.38263596012643665350944670744040
        a12=a10;    a13=a9;    a14=a8;    a15=a7;    a16=a6;    a17=a5
        a18=a4;    a19=a3;    a20=a2;    a21=a1
        gamma=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21]
            
    elseif r==10
         # 10th-order: (SS05)
         s=35
         a1=parse(BigFloat,"0.07879572252168641926390768")
         a2=parse(BigFloat,"0.31309610341510852776481247")
         a3=parse(BigFloat,"0.02791838323507806610952027")
         a4=parse(BigFloat,"-0.22959284159390709415121340")
         a5=parse(BigFloat,"0.13096206107716486317465686")
         a6=parse(BigFloat,"-0.26973340565451071434460973")
         a7=parse(BigFloat,"0.07497334315589143566613711")
         a8=parse(BigFloat,"0.11199342399981020488957508")
         a9=parse(BigFloat,"0.36613344954622675119314812")
         a10=parse(BigFloat,"-0.39910563013603589787862981")
         a11=parse(BigFloat,"0.10308739852747107731580277")
         a12=parse(BigFloat,"0.41143087395589023782070412")
         a13=parse(BigFloat,"-0.00486636058313526176219566")
         a14=parse(BigFloat,"-0.39203335370863990644808194")
         a15=parse(BigFloat,"0.05194250296244964703718290")
         a16=parse(BigFloat,"0.05066509075992449633587434")
         a17=parse(BigFloat,"0.04967437063972987905456880")
         a18=1-2*(a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15+a16+a17)
         a19=a17;    a20=a16;    a21=a15;    a22=a14;    a23=a13;    a24=a12
         a25=a11;    a26=a10;    a27=a9;    a28=a8;    a29=a7;    a30=a6
         a31=a5;    a32=a4;    a33=a3;    a34=a2;    a35=a1
         gamma=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18
             ,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35]

    else
        println("no method r")    
        @warn("There is not defined method of order r=$r")   
        s=0

    end

    if s!=0

        a=Array{Typeh}(undef,s+1)
        b=Array{Typeh}(undef,s)
        c=Array{Typeh}(undef,s)
        
        [b[i]=gamma[i]/2 for i in 1:s]
        [c[i]=gamma[i] for i in 1:s]
        a[1]=gamma[1]/2
        [a[i]=(gamma[i]+gamma[i-1])/2 for i in 2:s]
        a[s+1]=gamma[1]/2

    else

        a=[0]
        b=[0]
        c=[0]
    
    end


    return (s,a,b,c)


end


#
#  RKN splitting integrators
#     BM02: r=6, s=14 (14-stage symmetric 6th-order ABA)
#     BCE22: r=8, s=19 (19-stage symmetric 8th-order ABA)
#

function Splitting_rkn_Coefficients(r,h)

    Typeh=typeof(h)
    s=1

    if r==6
       # BM02 
       s=14

       a1=parse(BigFloat,"3.785931984061157e-2")
       a2=parse(BigFloat,"0.102635633102435")
       a3=parse(BigFloat,"-2.586788826655869e-2")
       a4=parse(BigFloat,"0.314241403071447")
       a5=parse(BigFloat,"-0.130144459517415")
       a6=parse(BigFloat,"0.106417700369543")
       a7=parse(BigFloat,"-8.794243128510581e-3")
       #a8=0.207305069056895
       a8=1-2*(a1+a2+a3+a4+a5+a6+a7)
   
       b1=parse(BigFloat,"9.171915262446165e-2")
       b2=parse(BigFloat,"0.183983170005006")
       b3=parse(BigFloat,"-5.653436583288827e-2")
       b4=parse(BigFloat,"4.914688774712854e-3")
       b5=parse(BigFloat,"0.143761127168358")
       b6=parse(BigFloat,"0.328567693746804")
       #b7=-0.196411466486454 
       b7=1/2-(b1+b2+b3+b4+b5+b6)
   
       a=[a1, a2, a3, a4, a5, a6, a7, a8, a7, a6, a5, a4, a3, a2, a1]
       c=[b1, b2, b3, b4, b5, b6, b7, b7, b6, b5, b4, b3, b2, b1,  0]
       b=zero(c)

    elseif r==8
        # BCE22
        s=19

        a1 = parse(BigFloat,"0.0505805") 
        b1 = parse(BigFloat,"0.129478606560536730662493794395")
        a2 = parse(BigFloat,"0.149999 ")
        b2 = parse(BigFloat,"0.222257260092671143423043559581")
        a3 = parse(BigFloat,"-0.0551795510771615573511026950361") 
        b3 = parse(BigFloat,"-0.0577514893325147204757023246320")
        a4 = parse(BigFloat,"0.423755898835337951482264998051") 
        b4 = parse(BigFloat,"-0.0578312262103924910221345032763")
        a5 = parse(BigFloat,"-0.213495353584659048059672194633") 
        b5 = parse(BigFloat,"0.103087297437175356747933252265")
        a6 = parse(BigFloat,"-0.0680769774574032619111630736274") 
        b6 = parse(BigFloat,"-0.140819612554090768205554103887")
        a7 = parse(BigFloat,"0.227917056974013435948887201671") 
        b7 = parse(BigFloat,"0.0234462603492826276699713718626")
        a8 = parse(BigFloat,"-0.235373619381058906524740047732") 
        b8 = parse(BigFloat,"0.134854517356684096617882205068")
        a9 = parse(BigFloat,"0.387413869179878047816794031058") 
        b9 = parse(BigFloat,"0.0287973821073779306345172160211") 
    
        a10=1/2-(a1+a2+a3+a4+a5+a6+a7+a8+a9)
        b10=1-2*(b1+b2+b3+b4+b5+b6+b7+b8+b9)
    
        a=[a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a10, a9, a8, a7, a6, a5, a4, a3, a2, a1]
        c=[b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b9, b8, b7, b6, b5, b4, b3, b2, b1,  0]
        b=zero(c)


    else
        println("no method r")    
        @warn("There is not defined method of order r=$r")   
        s=0
        a=[0]
        b=[0]
        c=[0]

    end


    return (s,convert.(Typeh,a),convert.(Typeh,b),convert.(Typeh,c))


end

