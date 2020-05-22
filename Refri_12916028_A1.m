%MODEL HIDRODINAMIKA 1-D TOPOGRAFI KONSTAN (MODUL III)-A1
%MUHAMMAD REFRI ANSYARI (12916028)

%% INPUT
%Inisialisasi Variabel

clear all;

L         = 2000;   %Panjang Kanal (m)
dx        = 100;    %Besar Grid Horizontal (m)
d         = 10;     %Topografi (Kedalaman) Kanal (m)
A         = 1.5;    %Amplitudo Gelombang Max (m)
T         = 12;     %Periode Gelombang (sekon)
g         = 10;     %Percepatan Gravitasi (m/s^2)
t         = 172800; %Lama Simulasi (sekon) - Total lama simulasi adalah 2 hari

%% Inisialisasi Parameter Gelombang 
C0        = sqrt(g*d);       % Kecepatan gelombang (perairan dangkal)
sigma     = 2*3.14/T;        % Frekuensi Sudut
k         = 2*3.14/(T*C0);   % Bilangan Gelombang

%% Syarat Kestabilan Courant-Freiderichs-Lewy (CFL)
dt        = dx/C0; 
%% Penentuan Jumlah Grid Waktu dan Ruang
imax=(L/dx); %Ruang
nmax=(t/dt); %Waktu
%% Nilai Awal dan Syarat Batas
% Nilai Awal
for i=1:imax
        %Nilai elevasi pada seluruh grid saat t=1
         elev(1,i)=A*cosd(k*i*dx);
        %Nilai kedalaman total pada seluruh grid saat t=1
         h(1,i)=d+elev(1,i);
        %Nilai kecepatan awal pada seluruh grid saat t=1
         u(1,i)=(A*C0/h(1,i))*cosd(k*((i*dx)+(dx/2)));
end
% Syarat Batas 
for n=2:nmax
         %Nilai elevasi pada grid batas di seluruh waktu
         elev(n,1)=A*cosd(sigma*n*dt);
         %Nilai kedalaman total pada grid batas di seluruh waktu
         h(n,1)=d+elev(n,1);
end
%% PROSES
        % Menghitung nilai u
for n=1:nmax-1
    for i=1:imax-1
        u(n+1,i)=u(n,i)-((g*(dt/2)/dx)*(elev(n,i+1)-elev(n,i)));
    end

        % Menghitung nilai elevasi dan kedalaman 
for i=2:imax
        elev(n+1,i)=elev(n,i)-((h(n,i)*(dt/2)/dx)*(u(n+1,i)-u(n+1,i-1)));
        h(n+1,i)=d+elev(n+1,i);
end

        % Penentuan Syarat Batas di Hilir (H baru) 
        u(n+1,imax)=(A*C0/(h(n+1,imax)))*cosd((k*L)-(sigma*(n+1)*dt));
end

%% OUTPUT
%% Plot Besar u (Sepanjang Waktu)
figure (1)
subplot(2,1,1);
plot (u(:,3),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 300 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (u(:,8),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 800 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

figure (2)
subplot(2,1,1);
plot (u(:,13),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 1300 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (u(:,19),'m');
title ('Kecepatan Aliran Hilir-Hulu Pada Grid 1900 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

%% Plot Besar u (Sepanjang Ruang)
figure (3) 
subplot(2,1,1);
plot (u(2160,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-6 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');


subplot(2,1,2);
plot (u(4680,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-13 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

figure(4)
subplot(2,1,1);
plot (u(7920,:),'c');
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-22 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');


subplot(2,1,2);
plot (u(16560,:));
title ('Kecepatan Aliran Hilir-Hulu Pada Jam ke-46 (Sepanjang Ruang)','fontweight','b');
ylabel ('u(m/s^2)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');


%% Plot Nilai Elevasi Kanal (Sepanjang waktu)
figure (5)
subplot(2,1,1);
plot (elev(:,3),'g');
ylim ([-4 4]);
xlim ([0 17280]);title ('Elevasi Hilir-Hulu Pada Grid 300 m (Sepanjang Waktu)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(:,8),'g');
title ('Elevasi Hilir-Hulu Pada Grid 800 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

figure (6)
subplot(2,1,1);
plot (elev(:,13),'g');
title ('Elevasi Hilir-Hulu Pada Grid 1300 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(:,19),'g');
title ('Elevasi Hilir-Hulu Pada Grid 1900 m (Sepanjang Waktu)','fontweight','b');
ylim ([-4 4]);
xlim ([0 17280]);
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

%% Plot Nilai Elevasi Kanal (Sepanjang Ruang)
figure (7)
subplot(2,1,1);
plot (elev(2160,:),'r');
title ('Elevasi Hilir-Hulu Pada Jam ke-6 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(4680,:),'r');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-13 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

figure(8)
subplot(2,1,1);
plot (elev(7920,:),'r');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-22 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (elev(16560,:),'r');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-46 (Sepanjang Ruang)','fontweight','b');
ylabel ('Elevasi (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

%% Plot Nilai H (Kedalaman) Sepanjang Waktu
figure (9)
subplot(2,1,1);
plot (h(:,3),'g');
ylim ([6 14]);
xlim ([0 17280]);title ('Kedalaman Hilir-Hulu Pada Grid 300 m (Sepanjang Waktu)','fontweight','b');
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (h(:,8),'g');
title ('Kedalaman Hilir-Hulu Pada Grid 800 m (Sepanjang Waktu)','fontweight','b');
ylim ([6 14]);
xlim ([0 17280]);
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

figure(10)
subplot(2,1,1);
plot (h(:,13),'g');
title ('Kedalaman Aliran Hilir-Hulu Pada Grid 1300 m (Sepanjang Waktu)','fontweight','b');
ylim ([6 14]);
xlim ([0 17280]);
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (h(:,19),'g');
title ('Kedalaman Aliran Hilir-Hulu Pada Grid 1900 m (Sepanjang Waktu)','fontweight','b');
ylim ([6 14]);
xlim ([0 17280]);
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Waktu (detik) x10','FontSize',8,'fontweight','b');

%% Plot Nilai H (Kedalaman) Sepanjang Ruang
figure (11)
subplot(2,1,1);
plot (h(2160,:),'y');
title ('Kedalaman Hilir-Hulu Pada Jam ke-6 (Sepanjang Ruang)','fontweight','b');
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (h(4680,:),'y');
title ('Elevasi Aliran Hilir-Hulu Pada Jam ke-13 (Sepanjang Ruang)','fontweight','b');
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

figure(12)
subplot(2,1,1);
plot (h(7920,:),'y');
title ('Kedalaman Aliran Hilir-Hulu Pada Jam ke-22 (Sepanjang Ruang)','fontweight','b');
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');

subplot(2,1,2);
plot (h(16560,:),'y');
title ('Kedalaman Aliran Hilir-Hulu Pada Jam ke-46 (Sepanjang Ruang)','fontweight','b');
ylabel ('Kedalaman (m)','FontSize',12,'fontweight','b');
xlabel ('Panjang Kanal (m) x100','FontSize',8,'fontweight','b');