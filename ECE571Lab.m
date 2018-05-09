function main()
lab8()
lab9()
lab10()
lab11()
lab12()
end
function lab8()
v38 = [520e-6,2.34e-3,10.36e-3,20.26e-3,39.97e-3,49.6e-3];
i45 = [1e-6,10e-6,50e-6,100e-6,.2e-3,.25e-3];
v56 = [88.7e-3,.316,.626,.82,1.05,1.15,1.29];
i37 = [1e-6,10e-6,50e-6,100e-6,200e-6,.25e-3,.35e-3];
L = 3300e-6;
v12 = [.15e-3,.245e-3,.643e-3,1.13e-3,2.12e-3,3.1e-3,3.6e-3];
i12 = [1e-6,10e-6,50e-6,100e-6,.2e-3,.3e-3,.35e-3];
v312 = [212e-3,1.23,2.05,2.46,2.92,3.2,3.4,3.48];
i312 = [1e-6,10e-6,50e-6,100e-6,.2e-3,.3e-3,.4e-3,.45e-3];
v110 = [.948,1.9,2.8,3.35,3.8,5.02,4.48,4.67];
i110 = [1e-6,10e-6,50e-6,100e-6,.2e-3,.3e-3,.4e-3,.45e-3];
rsfit = polyfit(i45,v38,1)
rs = (pi / log(2)) * rsfit(1)
rwfit = polyfit(i37,v56,1)
w = rs * L / rwfit(1)
%r12 = v12 ./ i12
%r312 = v312 ./ i312
%r110 = v110 ./ i110
r12 = polyfit(i12,v12,1)
r312 = polyfit(i312,v312,1)
r110 = polyfit(i110,v110,1)
rAl = r12(1) / 2
rB = r312(1) / 4
rC = (r110(1) - r12(1) - 2 * rB) / 2
%figure
%scatter(i12,v12)
%figure
%scatter(i312,v312)
%figure
%scatter(i110,v110)
%figure
%scatter(i45,v38)
end
function lab9()
sSmall = .3e-3;
sBig = .8e-3;
ivV = linspace(.1,1.9,19);
%subtract off leakage current
smallLeak = -16.5e-9;
bigLeak = -90.5e-9;
currentSmallDiode = [.14e-6,1.64e-6,14.6e-6,79.6e-6,.195e-3,.317e-3,.421e-3,.509e-3,.579e-3,.638e-3,.692e-3,.747e-3,.805e-3,.865e-3,.929e-3,.995e-3,1.065e-3,1.138e-3,1.213e-3];
currentDensitySmall = (currentSmallDiode - smallLeak) ./ (sSmall * sSmall);
jsSmall = polyfit(ivV(1:4),log(currentDensitySmall(1:4)),1)
currentBigDiode = [4.5e-9,1.73e-6,32.12e-6,.1484e-3,.289e-3,.411e-3,.507e-3,.581e-3,.642e-3,.698e-3,.755e-3,.817e-3,.881e-3,.948e-3,1.02e-3,1.095e-3,1.173e-3,1.255e-3,1.340e-3];
currentDensityBig = (currentBigDiode - bigLeak) ./ (sBig * sBig);
jsBig = polyfit(ivV(1:4),log(currentDensityBig(1:4)),1)
figure
jSmall = semilogy(ivV,currentDensitySmall);
jSmall.Marker = 'o';
hold on
%ySmall =  exp(jsSmall(1) .* ivV(1:4) + jsSmall(2));
idealSmall = semilogy(ivV(1:5),logPlot(ivV(1:5),jsSmall(1),jsSmall(2)));
hold off
figure
jBig = semilogy(ivV,currentDensityBig);
jBig.Marker = 'o';
hold on
%yBig =  exp(jsBig(1) .* ivV(1:4) + jsBig(2));
idealBig = semilogy(ivV(1:5),logPlot(ivV(1:5),jsBig(1),jsBig(2)));
hold off
cvV = linspace(-20,-2,10);
cSmallDiode1k = 10 ^ -12 .* [8.15058,8.27375,8.64152,9.09706,9.78703,12.1729,16.0637,17.7857,22.0533,33.7795];
cBigDiode1k = 10 ^ -12 .* [8.41612,16.6969,20.9739,20.9626,19.9782,20.9313,23.2835,27.8592,38.5609,72.7088];
cSmallDiode100k = 10 ^ -12 .* [7.7757,8.1448,8.62668,9.20643,9.93123,10.6531,11.6322,13.7460,19.1400,32.4739];
cBigDiode100k = 10 ^ -12 .* [14.5985,15.953,16.6940,17.7655,17.4618,20.4370,23.0042,27.4131,37.4822,69.1313];
sSD1k = plot1c2(cvV,cSmallDiode1k)
sSD100k = plot1c2(cvV,cSmallDiode100k)
sBD1k = plot1c2(cvV,cBigDiode1k)
sBD100k = plot1c2(cvV,cBigDiode100k)
ndS1k = ndify(sSmall,sSD1k(1))
ndS100k = ndify(sSmall,sSD100k(1))
ndB1k = ndify(sBig,sBD1k(1))
ndB100k = ndify(sBig,sBD100k(1))
vBiS1k = findVbi(ndS1k,sSD1k(2),sSmall)
vBiS100k = findVbi(ndS100k,sSD100k(2),sSmall)
vBiB1k = findVbi(ndB1k,sBD1k(2),sBig)
vBiB100k = findVbi(ndB100k,sBD100k(2),sBig)
nd = 1e15;
nSmall = findN(jsSmall(1))
nBig = findN(jsBig(1))
lpSmall = lp(exp(jsSmall(2)),nd)
lpBig = lp(exp(jsBig(2)),nd)
tsS = 1.48e-9;
tsB = 1.88e-9;
vFS = 3.8;
vRS = 4.7;
vFB = 3.14;
vRB = 3.18;
tpS = findtp(tsS,vFS,vRS)
tpB = findtp(tsB,vFB,vRB)
lpTS = lpTrans(tpS)
lpTB = lpTrans(tpB)
end
function y = findtp(ts,vF,vR)
y = ts * log(1 + (vF / vR));
end
function y = lpTrans(tp)
dp = 10;
y = sqrt(dp * tp);
end
function y = ndify(s,cs)
q = 1.6021766208e-19;
eSi = 11.68 * 8.854187817e-12;
y = -1 / ( .5 * q * cs * eSi * ((s ^ 2) ^ 2));
end
function y = findVbi(nd,ci,s)
q = 1.6021766208e-19;
eSi = 11.68 * 8.854187817e-12;
y = ci * .5 * q * nd * eSi * ((s ^ 2) ^ 2);
end
function y = findN(jsExp)
y = 1 / (jsExp * .026);
end
function y = plot1c2(v,c)
figure
c2 = 1 ./ (c .^ 2);
p = plot(v,c2);
p.Marker = 'o';
y = polyfit(v,c2,1);
end
function y = logPlot(x,m,b)
y =  exp(m * x + b);
end
function y = lp(js,nd)
dp = 10;
ni2 = 1.5e10 ^ 2;
q = 1.6021766208e-19;
%nd = 5;
y = q * dp * ni2 / (js * nd);
end
function lab10()
sSmall = .5e-3;
sBig = 1e-3;
v = linspace(-5,4,10);
tox = [];
cox = [];
cinv = [];
coxp = [];
cinvp = [];
% gS1k;gS100k;fS1k;fS100k;gB1k;gB100k;fB1k;fB100k;
c = 10 ^ -12 .* [117.0518,116.4584,111.7531,105.6328,97.9698,107.7867,113.1229,115.9079,117.1398,118.4676;21.2072,21.6059,22.4896,24.8516,51.7554,92.2331,103.8094,114.3687,116.3698,117.4887;37.4086,35.72,32.1565,34.1872,34.6668,38.8543,39.0118,38.0741,37.3332,35.4458;12.1916,12.3088,14.1208,25.5074,36.2222,35.5443,38.5405,36.4878,36.6084,36.9453;379.588,378.823,374.628,333.156,294.552,321.225,370.822,376.739,380.838,381.159;47.3857,46.613,47.0898,60.7291,171.665,296.839,366.402,372.109,374.471,377.681;140.615,137.728,122.44,118.641,129.393,135.881,138.824,141.059,140.4,141.55;42.4329,40.3633,56.339,97.49,127.92,138.32,139.018,138.584,140.225,141.529];
for i = 1:4
    tox = [tox;toxify(sSmall,c(i,10))];
    coxp = [coxp;c(i,10) / (sSmall ^ 2)];
    cox = [cox;c(i,10)];
    if(mod(i,2) == 0)
%        cinv = [cinv;c(i,1) / (sSmall ^ 2)];
%        cinv = [cinv;c(i,1) / (sSmall ^ 2)];
        cinv = [cinv;c(i,1)];
        cinv = [cinv;c(i,1)];
    end
end
for i = 5:8
    tox = [tox;toxify(sBig,c(i,10))];
    coxp = [coxp;c(i,10) / (sBig ^ 2)];
    cox = [cox;c(i,10)];
    if(mod(i,2) == 0)
%        cinv = [cinv;c(i,1) / (sBig ^ 2)];
%        cinv = [cinv;c(i,1) / (sBig ^ 2)];
        cinv = [cinv;c(i,1)];
        cinv = [cinv;c(i,1)];
    end
end
display(tox)
display(cox)
display(cinv)
nD = 10 ^ 16 .* 10 ^ 6 .* ones(8,1);
nDOld = ones(8,1);
phiOld = ones(8,1);
phi = phiOld;
cdep = 1 ./ ((1 ./ cinv) - (1 ./ cox))
index = 0;
while ~isequal(nDOld,nD) || ~isequal(phiOld,phi)
    nDOld = nD;
    phiOld = phi;
    phi = abs(calcPhi(nDOld));
    nD = abs(calcND(phi,cdep));
    index = index + 1;
    if index > 10000
        break
    end
end
%display(nD - nDOld)
display(nD)
%display(phi - phiOld)
display(phi)
for index = 0:3
    plotCV(v,c,1 + 2 * index)
end
cdepFB = calcCDep(nD)
cFBp = coxp .* cdepFB ./ (coxp + cdepFB)
cFB = [];
for index = 1:4
    cFB = [cFB;undop(sSmall,cFBp(index,:))];
end
for index = 5:8
    cFB = [cFB;undop(sBig,cFBp(index,:))];
end
display(cFB)
end
function y = undop(s,cp)
y = (s ^ 2) .* cp;
end
function y = calcCDep(nD)
q = 1.6021766208e-19;
eox = 3.9;
e0 = 8.854187817e-12;
eSi = eox * e0;
k = 1.38064852e-23;
T = 300;
y = sqrt(nD .* (q ^ 2 * eSi / (k * T)));
end
function plotCV(v,c,index)
figure
plot(v,c(index,:),'o-');
hold on
plot(v,c(index + 1,:),'o-')
hold off
end
function y = calcND(phi,cDep)
q = 1.6021766208e-19;
eox = 3.9;
e0 = 8.854187817e-12;
eSi = eox * e0;
y = phi .* 4 .* cDep .^ 2 ./ (q .* eSi);
end
function y = calcPhi(nD)
kTq = .02586;
ni = 1.5e10 * 10 ^ 6;
y = (2 .* kTq) .* log(nD ./ ni);
end
function y = toxify(s,cAcc)
a = s ^ 2;
eOx = 3.9 * 8.854187817e-12;
y = eOx * a ./ cAcc;
end
function lab11()
iD2fit = [];
vDS = -16;
vGS = linspace(-1,-5,5);
l = 10 ^ -6 .* [25;50;100;200];
iD = 10 ^ -6 .* [-1.9082,-3.7656,-273,-862,-1912;-3.341,-59.568,-279.9,-729.02,-1131.6;-2.7278,-24.327,-155.63,-322.50,-543.8;-1.6367,-11.342,-60.998,-180.96,-286.79];
for i = 1:4
    iD2fit = [iD2fit;plotiD2(vGS,iD(i,:))];
end
vT = iD2fit(:,2) ./ iD2fit(:,1)
iD1 = 1 ./ iD;
lRealfit = [];
dL = [];
figure
hold on
for i = 1:length(iD1)
    p = plot(l(:,1),iD1(:,i));
    p.Marker = 'o';
    lRealfit = [lRealfit;polyfit(l(:,1),iD1(:,i),1)];
    plot(l(:,1),lRealfit(i,1) .* l(:,1) + lRealfit(i,2));
    dL = [dL;lRealfit(i,2) ./ lRealfit(i,1)];
end
hold off
display(lRealfit)
display(dL)
meanDL = mean(dL)
standardDevDL = std(dL)
lReal = l - meanDL
eO = 8.854187817e-12;
eSi = 3.9;
tOx = .68e-6
cOx = eO * eSi / tOx;
w = 1000e-6;
upSi = 1500;
up = 10 ^ 4 .* (iD2fit(:,1) .^ 2) .* 2 .* l ./ (cOx .* w)
meanup = mean(up)
upPercent = 100 - (100 * meanup / upSi)
end
function y = plotiD2(v,iD)
figure
iD2 = sqrt(abs(iD));
p = plot(v,iD2);
p.Marker = 'o';
y = polyfit(v,iD2,1);
hold on
plot(v,y(1) .* v + y(2));
hold off
end
function lab12()
rs = 15;
r = 100e3;
s = 25e-6;
rtLength = r * s / rs
invNum = 10 ^ 6;
aInvNum = rtLength * s * invNum
v = 10;
up = .15;
power = v ^ 2 / r * invNum
vsup = linspace(-11,-15,5);
lS = 10 ^ -6 .* [25;50;100;200];
wS = 10 ^ -6 .* [125;250;500;1000];
lL = 10 ^ -6 .* [125;250;500;1000];
wL = 10 ^ -6 .* [25;50;100;200];
current = -10 ^ -3 .* [.2248,.3547,.4206,.5301,.6803;.1491,.23743,.2669,.3494,.3954;.2055,.2665,.3101,.3768,.4299;.2248,.2743,.3149,.3442,.4146];
period = 10 ^ -6 .* [9.92,8.60,7.78,6.47,4.11;26.67,20.48,17.77,13.66,10.58;44.65,39.13,36.25,31.71,27.50;111.05,92.89,82.30,79.98,69.17];
frequency = 1 ./ period
freqTheory = up .* abs(vsup) ./ (350 .* lS .^ 2)
propDelay = period ./ 10
ringPower = 5 .* vsup .* current
%delayVPowerfit = [];
figure
for i = 1:4
    loglog(ringPower(i,:),propDelay(i,:),'o-');
 %   delayVPowerfit = [delayVPowerfit;]
    loglog(ringPower(i,:),i .* 10 ^ -6 .* exp(-1 .* ringPower(i,:)))
    hold on
end
hold off
end