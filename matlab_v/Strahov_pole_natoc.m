clc, clear, format longg
close all

% mereni pole 16.09.2019
bod1 = [50.0792144685 14.385311273]; %pocatek (jiznejsi bod)
bod2 = [50.0794838003 14.385289759]; %urceni kroku

[vzdal, aaa, aaaa, aaaaa] = gedistance(bod1(1), bod1(2), bod2(1), bod2(2))

pocet_bodu = 750000;

diffLAT = bod2(1) - bod1(1);
diffLON = bod2(2) - bod1(2);

deltaLAT = diffLAT / vzdal;
deltaLON = diffLON / vzdal;

azimut = azimuth(bod1(1), bod1(2), bod2(1), bod2(2))
alfa = 360 - azimut;

R1 = 6378000; ;
R2 = cosd(bod1(1, 1)) * R1;

deltalon = (sind(azimut - 90) / R2) * (180 / pi);
deltalat = (cosd(azimut - 90) / R1) * (180 / pi);

%%======================================

%{

bodyOSA=[bod1, 10 ;zeros(pocet_bodu,3)];

for i=2:(pocet_bodu+1)

    bodyOSA(i,1)=bodyOSA(i-1,1)+deltaLAT;
    bodyOSA(i,2)=bodyOSA(i-1,2)+deltaLON;
    bodyOSA(i,3)=10 + (i-1)*5;

end

bodyP1=[bodyOSA(:,1) -   deltalat,bodyOSA(:,2) -   deltalon,bodyOSA(:,3)-5];
bodyP2=[bodyOSA(:,1) - 2*deltalat,bodyOSA(:,2) - 2*deltalon,bodyOSA(:,3)-10];
bodyP3=[bodyOSA(:,1) - 3*deltalat,bodyOSA(:,2) - 3*deltalon,bodyOSA(:,3)-15];
bodyL1=[bodyOSA(:,1) +   deltalat,bodyOSA(:,2) +   deltalon,bodyOSA(:,3)+5];
bodyL2=[bodyOSA(:,1) + 2*deltalat,bodyOSA(:,2) + 2*deltalon,bodyOSA(:,3)+10];
bodyL3=[bodyOSA(:,1) + 3*deltalat,bodyOSA(:,2) + 3*deltalon,bodyOSA(:,3)+15];

ss=[bodyOSA;bodyL1;bodyL2;bodyL3;bodyP1;bodyP2;bodyP3];
ss_model_rovina = ss;
ss_model_rovina(:,3) = 400;
ss_model_sklon = ss;
ss_model_sklon(:,3) = 400 - ss(:,3)/100;% 5% spád

bodyOsaVedle=[bodyOSA(1:end,1)-4*deltalat,bodyOSA(1:end,2)-4*deltalon]; % 4 m vlevo

ss_reality_rovina = ss;
ss_reality_rovina(:,3) = 400 + ss(:,3)/1000;
ss_reality_sklon = ss;
ss_reality_sklon(:,3) = ss_model_sklon(:,3) + ss(:,3)/1000;

disp('dopoctene body')

% polygon

plg = [bodyP3(1,1)   bodyP3(1,2)
       bodyP3(end,1) bodyP3(end,2)
       bodyL3(end,1) bodyL3(end,2)
       bodyL3(1,1)   bodyL3(1,2)
       ];

disp('transformace do UTM')

ss_UTM = ss;
bodyOSA_UTM = bodyOSA;
bodyOsaVedle_UTM = bodyOsaVedle;
ss_model_rovina_UTM = ss_model_rovina;
ss_model_sklon_UTM = ss_model_sklon;
ss_reality_sklon_UTM = ss_reality_sklon;
ss_reality_rovina_UTM = ss_reality_rovina;
plg_UTM = plg;

[ss_UTM(:,1),ss_UTM(:,2)] = ll2utm(ss(:,1:2));
[bodyOSA_UTM(:,1),bodyOSA_UTM(:,2)] = ll2utm(bodyOSA(:,1:2));
[bodyOsaVedle_UTM(:,1),bodyOsaVedle_UTM(:,2)] = ll2utm(bodyOsaVedle(:,1:2));
[ss_model_rovina_UTM(:,1),ss_model_rovina_UTM(:,2)] = ll2utm(ss_model_rovina(:,1:2));
[ss_model_sklon_UTM(:,1),ss_model_sklon_UTM(:,2)] = ll2utm(ss_model_sklon(:,1:2));
[ss_reality_sklon_UTM(:,1),ss_reality_sklon_UTM(:,2)] = ll2utm(ss_reality_sklon(:,1:2));
[ss_reality_rovina_UTM(:,1),ss_reality_rovina_UTM(:,2)] = ll2utm(ss_reality_rovina(:,1:2));
[plg_UTM(:,1),plg_UTM(:,2)] = ll2utm(plg(:,1:2));

disp('zapis do souboruuu')

dlmwrite( ['strah_dtm.txt'], ss_UTM, 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_osa.txt'], bodyOSA_UTM(:,1:2), 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_osaVedle.txt'], bodyOsaVedle_UTM, 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_design_rovina.txt'], ss_model_rovina_UTM, 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_design_sklon.txt'], ss_model_sklon_UTM, 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_reality_sklon.txt'], ss_reality_sklon_UTM, 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_reality_rovina.txt'], ss_reality_rovina_UTM, 'delimiter', ',', 'precision', '%.9f' );
dlmwrite( ['strah_vnejsi_plg.txt'], plg_UTM, 'delimiter', ',', 'precision', '%.9f' );

disp('zapsano do souboruuu')

%}

%{
 %% kontrola vzdalenosti
%vzd = gedistance(0,0,deltalat,deltalon)
vzdd =  gedistance(bodyOSA(1,1),bodyOSA(1,2),bodyL1(1,1),bodyL1(1,2));
vzd1 =  gedistance(bodyOSA(1,1),bodyOSA(1,2),bodyOSA(2,1),bodyOSA(2,2));
vzd2 =  gedistance(bodyL2(1,1),bodyL2(1,2),bodyOSA(1,1),bodyOSA(1,2));
vzd3 =  gedistance(bodyL1(1,1),bodyL1(1,2),bodyL1(2,1),bodyL1(2,2));
vzd4 =  gedistance(bodyL1(1,1),bodyL1(1,2),bodyL2(1,1),bodyL2(1,2));

%ss1=[
%50.079156097 14.385292642
%50.079174697 14.385292675
%];

%
plot(ss(:,2),ss(:,1),'r*')
hold on
axis equal
plot(bodyOSA(:,2),bodyOSA(:,1),'b*')
plot(bodyP2(:,2),bodyP2(:,1),'k*')
%plot(ss(1,2),ss(1,1),'g*')
plot(plg(:,2),plg(:,1),'g*')

%plot(bodyL2(1,2),bodyL2(1,1),'g*')
%plot(bodyL1(1,2),bodyL1(1,1),'g*')
%plot(bodyOsaVedle(:,2),bodyOsaVedle(:,1),'g*')

%}
