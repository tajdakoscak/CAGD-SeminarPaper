function D = bezier(B,t)
% Opis:
% bezier vrne točke na Bezierjevi krivulji pri danih parametrih
%
% Definicija:
% b = bezier(B,t)
%
% Vhodna podatka:
% B matrika velikosti n+1 x d, ki predstavlja kontrolne točke
% Bezierjeve krivulje stopnje n v d-dimenzionalnem prostoru,
% t seznam parametrov dolžine k, pri katerih računamo vrednost
% Bezierjeve krivulje
%
% Izhodni podatek:
% b matrika velikosti k x d, kjer i-ta vrstica predstavlja točko
% na Bezierjevi krivulji pri parametru iz t na i-tem mestu

k = length(t);
[~,d] = size(B);
D = [];

for i = 1:k % tok točk moramo zračunat
    tocka = [];
    for j=1:d % dimenzija točk 
        shema = decasteljau(B(:,j), t(i));
        [st_vrstic, st_stolpcev] = size(shema);
        if st_vrstic > 0
            koordinata = shema(1,st_stolpcev); 
            tocka = [tocka, koordinata];
        end
    end
    D = [D; tocka];

end

end
