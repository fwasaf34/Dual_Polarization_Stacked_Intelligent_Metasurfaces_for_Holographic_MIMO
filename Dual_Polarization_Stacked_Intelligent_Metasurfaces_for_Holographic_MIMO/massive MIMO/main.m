clc
clear
close all


%%
count=0;
for PdBm=[-10:5:40]
    count=count+1;
    for seed =1:1000
        rng(seed)

        %%EM
        c=3*10^8;
        f=28*10^9;
        lambda=c/f;

        %%TX
        zTx=0;

        %%RX
        zRx=250;


        %%noise
        noisedBm=-110;
        noise=db2pow(noisedBm)*10^-3;

        %%P
        P=db2pow(PdBm)*10^-3;

        %%algorithm
        alp=0.4;

        %%channel
        rx=32;
        tx=256;

        GTilde=zeros(rx, tx);

        PL=20*log10(4*pi/lambda)+10*3.5*log10(zRx)+9;
        rho_sq=1/db2pow(PL);

        GTilde = sqrt((1-alp)*rho_sq/2) * (randn(rx,tx) + 1j*randn(rx,tx));

        G=GTilde;

        SE(seed) = calculate_spectral_efficiency(G, P, noise);
        EE(seed) = calculate_spectral_efficiency(G, P, noise)/P;

    end
    SEme(count)=mean(SE)
    EEme(count)=mean(EE)
end



