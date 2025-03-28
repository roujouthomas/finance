%% Script complet pour le pricing d'une Put Américaine
% Méthode Crank-Nicolson + PSOR
clear; close all; clc;

% Définition des paramètres de l'option (exemple)
r = 0.06;        % taux d'intérêt
sigma = 0.3;     % volatilité
T = 1;           % échéance (en années)
K = 10;          % strike de l'option
Smax = 15;       % valeur maximale de S (grille spatiale)
NS = 100;        % nombre de sous-intervalles spatiaux
Nt = 800;        % nombre de pas temporels
omega = 1.2;     % paramètre de relaxation PSOR

% Appel de la fonction de pricing
[S, V] = BS_Crank_Nicolsons_PSOR(r, sigma, T, K, Smax, NS, Nt, omega);

%% Affichage complémentaire (optionnel)
% La fonction BS_Crank_Nicolsons_PSOR réalise déjà l'affichage.
% Vous pouvez ajouter ici d'autres analyses ou tracés si besoin.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction locale : Pricing d'une Put Américaine par Crank-Nicolson + PSOR
function [S, V] = BS_Crank_Nicolsons_PSOR(r, sigma, T, K, Smax, NS, Nt, omega)
    % Discrétisation spatiale et temporelle
    dS = Smax / NS;
    dt = T / Nt;
    S = linspace(0, Smax, NS+1)'; % grille spatiale (vecteur colonne)
    
    % Initialisation de la matrice solution V :
    % Chaque colonne correspond à un instant de temps, t = T est la condition terminale.
    V = zeros(NS+1, Nt+1);
    % Condition terminale (à t = T) : payoff d'une Put = max(K - S, 0)
    V(:, end) = max(K - S, 0);
    
    % Conditions aux limites pour tout t :
    V(1, :) = K;    % à S = 0, l'option vaut K (exercice optimal)
    V(end, :) = 0;  % à S = Smax, l'option vaut 0
    
    % Nombre de points intérieurs (excluant S = 0 et S = Smax)
    M = NS - 1;
    
    % Calcul des coefficients pour le schéma Crank-Nicolson
    % (cf. pages 157-163 de aulas_alunos.pdf)
    alpha = zeros(M,1); 
    beta  = zeros(M,1); 
    gamma = zeros(M,1);
    for i = 1:M
        % Ici, i correspond au point S(i+1) (puisque S(1)=0)
        alpha(i) = 0.25 * dt * (sigma^2 * (i)^2 - r * i);
        beta(i)  = -0.5 * dt * (sigma^2 * (i)^2 + r);
        gamma(i) = 0.25 * dt * (sigma^2 * (i)^2 + r * i);
    end
    
    % Construction de la matrice A (partie implicite) pour les points intérieurs (taille M x M)
    A = zeros(M, M);
    for i = 1:M
        A(i,i) = 1 - beta(i);
        if i > 1
            A(i, i-1) = -alpha(i);
        end
        if i < M
            A(i, i+1) = -gamma(i);
        end
    end
    
    % Construction de la matrice B (partie explicite) pour la partie droite
    B = zeros(M, M);
    for i = 1:M
        B(i,i) = 1 + beta(i);
        if i > 1
            B(i, i-1) = alpha(i);
        end
        if i < M
            B(i, i+1) = gamma(i);
        end
    end
    
    % Paramètres pour la méthode PSOR
    tol = 1e-6;
    max_iter = 10000;
    
    % Récurrence temporelle : on remonte de n = Nt à n = 1 (intégration en temps inverse)
    for n = Nt:-1:1
        % Construction du vecteur d pour la résolution à l'instant n
        d = B * V(2:NS, n+1);
        % Ajustement avec les conditions aux limites
        d(1) = d(1) + alpha(1) * V(1, n+1);   % à S = 0 : V = K
        d(end) = d(end) + gamma(end) * V(end, n+1); % à S = Smax : V = 0
        
        % Définition du payoff pour les points intérieurs
        payoff = max(K - S(2:NS), 0);
        
        % Initialisation de la solution pour PSOR par la valeur de l'instant n+1
        U_guess = V(2:NS, n+1);
        U_new = U_guess;
        
        error_PSOR = inf;
        iter = 0;
        while error_PSOR > tol && iter < max_iter
            U_old = U_new;
            for i = 1:M
                sum_term = 0;
                if i > 1
                    sum_term = sum_term + A(i, i-1) * U_new(i-1);
                end
                if i < M
                    sum_term = sum_term + A(i, i+1) * U_new(i+1);
                end
                U_temp = (d(i) - sum_term) / A(i,i);
                % Relaxation PSOR
                U_temp = U_new(i) + omega * (U_temp - U_new(i));
                % Projection pour respecter la contrainte (U >= payoff)
                U_new(i) = max(U_temp, payoff(i));
            end
            error_PSOR = norm(U_new - U_old, inf);
            iter = iter + 1;
        end
        
        % Stockage de la solution pour l'instant n
        V(2:NS, n) = U_new;
        V(1, n) = K;
        V(NS+1, n) = 0;
    end
    
    % Affichage de la solution à t = 0 (valeur de l'option initiale)
    figure;
    plot(S, V(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(S, max(K - S,0), 'r--', 'LineWidth', 2);
    xlabel('S'); ylabel('Prix de l''option');
    title('Option Put Américaine : Prix à t = 0 vs. Payoff');
    legend('Prix calculé', 'Payoff');
    grid on;
    
    % Affichage de la surface de la solution V(S,t)
    figure;
    tgrid = linspace(0, T, Nt+1);
    surf(S, tgrid, V', 'EdgeColor', 'none');
    xlabel('S'); ylabel('t'); zlabel('Prix de l''option');
    title('Surface de la solution pour la Put Américaine (Crank-Nicolson + PSOR)');
    colorbar;
    shading interp;
end
