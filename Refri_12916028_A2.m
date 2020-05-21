%MODEL HIDRODINAMIKA 1-D TOPOGRAFI PERSAMAAN TRANSPORT (MODUL III)-A2
%MUHAMMAD REFRI ANSYARI (12916028)

%% INPUT
%Inisialisasi Variabel
clear all;

L         = 2000;           %Panjang Kanal (m)
dx        = 100;            %Besar Grid Horizontal (m)
A         = 1.5;            %Amplitudo Gelombang Max (m)
T         = 12;             %Periode Gelombang (sekon)
g         = 10;             %Percepatan Gravitasi (m/s^2)
t         = 86400;          %Lama Simulasi (sekon) - Total lama simulasi adalah 1 hari
h         = 12:-7/20:3;     %Interval perubahan kedalaman ( hilir-hulu)
dhilir    = 12              %Kedalaman di Hilir
dhulu     = 3               %Kedalaman di Hulu
hx        = 1:50:2000;
%% Inisialisasi Parameter Gelombang
C0        = sqrt(g*dhilir); %Kecepatan gelombang (perairan dangkal)
sigma     = 2*3.14/T;       %Frekuensi Sudut
k         = 2*3.14/(T*C0);  %Bilangan Gelombang
%% Penentuan Langkah Waktu
dt        =5; 
%% Penentuan Jumlah Grid Waktu dan Ruang
imax=(L/dx); %Ruang
nmax=(t/dt); %Waktu
%% Nilai Awal dan Syarat Batas
%Nilai Awal
for i=1:imax
        % Nilai elevasi pada seluruh grid saat t=1
        elev(1,i)=A*cos(k*i*dx);
        %Nilai kecepatan awal pada seluruh grid saat t=1
        u(1,i)=A*C0*cos(k*(hx(i)+(0.5*dx)));
end
%Syarat Batas 
for n=1:nmax
        %Nilai elevasi pada grid batas di seluruh waktu
        elev(n,1)=A*cos(sigma*n*dt);
        %Nilai kedalaman total pada grid batas di seluruh waktu
        u(n,imax)=A*C0*cos(k*L-(sigma*n*dt));
end

        u(2,:)=u(1,:);
        elev(2,:)=elev(1,:);
%% PROSES
            % Menghitung nilai u
for n=2:nmax-1
    for i=2:imax-1
            u(n+1,i)=u(n-1,i)-(0.5*g*(dt/dx)*(h(i+1))+h(i)*(elev(n,i+1)-elev(n,i)));
    end
    
            % Menghitung nilai elevasi dan kedalaman 
    for i=2:imax-1
             elev(n+1,i)=elev(n,i)-(dt/dx)*(u(n+1,i)-u(n+1,i-1));     
    end
    
            % Penentuan Syarat Batas di Hilir (H baru)
    for i=2:imax-1
            u(n+1,i)=u(n,i)-(0.5*g*(dt/dx)*(h(i+1))+h(i)*(elev(n+1,i+1)-elev(n+1,i)));
    end
end

%%
% %Output Hasil
%% Plot Besar u (Sepanjang Waktu)
figure (1)
subplot(2,1,1);
plot (u(:,3),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 300 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (u(:,8),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 800 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

figure (2)
subplot(2,1,1);
plot (u(:,13),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 1300 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (u(:,19),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 1900 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

%% Plot Besar u (Sepanjang Ruang)
figure (3) 
subplot(2,1,1);
plot (u(4320,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-6 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (u(9360,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-13 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

figure(4)
subplot(2,1,1);
plot (u(12960,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-18 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (u(15840,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-22 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');


%% Plot Nilai Elevasi Kanal (Sepanjang waktu)
figure (5)
subplot(2,1,1);
plot (elev(:,3),'g');
title ('Elevasi Hilir-Hulu Pada Grid 300 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(:,8),'g');
title ('Elevasi Hilir-Hulu Pada Grid 800 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

figure (6)
subplot(2,1,1);
plot (elev(:,13),'g');
title ('Elevasi Hilir-Hulu Pada Grid 1300 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(:,19),'g');
title ('Elevasi Hilir-Hulu Pada Grid 1900 m (Sepanjang Waktu)','fontweight','b');
xlim ([0 17280])
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

%% Plot Nilai Elevasi Kanal (Sepanjang Ruang)
figure (7)
subplot(2,1,1);
plot (elev(4320,:),'r');
title ('Elevasi Hilir-Hulu Pada Jam ke-6 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(9360,:),'r');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-13 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

figure(8)
subplot(2,1,1);
plot (elev(12960,:),'r');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-18 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(15840,:),'r');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-22 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');
