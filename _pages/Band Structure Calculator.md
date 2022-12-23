---
title: "Band Structure Calculator based on MATLAB"
layout: collection
permalink: /projects/BandStructureCalculator/
entries_layout: grid
classes: wide
excerpt: "Hamiltonian Matrix based Kronig-Penny Model Simulation"
---

```
function E = KPhex_6( N, v, d, iter )
E = zeros((iter+1),(iter+1),(N+1)*(N+1));

% Converts {1,2,3,4,...} -> {0,1,-1,2,-2,...}
j = 1:N+1;
a = (1 + (2*j - 1).*(-1).^j)./4;

% Naive encoding, using meshgrid to generate pairwise combinations
[p,q] = meshgrid(a,a);
n = [p(:) q(:)];

% Increasing energy encoding
n2 = zeros((N+1)*(N+1),1);
n = [n2,n];

% Appending n^2 column to sort by
n(:,1) = n(:,2).*n(:,2) + n(:,3).*n(:,3);

% Sorting by n^2 to give increasing energy encoding
[T,I] = sort(n(:,1));
n = n(I,:);
nx = n(:,2); ny = n(:,3); n2 = n(:,1);
h = zeros((N+1)^2, (N+1)^2);

% Populating off-diagonal elements
for i = 1:(N+1)*(N+1)
    for j = i+1:(N+1)*(N+1)

        if nx(i) == nx(j)
            xpart1 = 2*d/3;
            xpart2 = 2*d/3;
            xpart3 = 2*d/3;
            xpart4 = 2*d/3;
        else
            m = nx(j) - nx(i);
            xpart1 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(1/6-d/3)) + exp(1i*2*pi*(m)*(1/6+d/3)));
            xpart2 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(2/6-d/3)) + exp(1i*2*pi*(m)*(2/6+d/3)));
            xpart3 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(4/6-d/3)) + exp(1i*2*pi*(m)*(4/6+d/3)));
            xpart4 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(5/6-d/3)) + exp(1i*2*pi*(m)*(5/6+d/3)));
        end

        if ny(i) == ny(j)
            ypart1 = 2*d/sqrt(3);
            ypart2 = 2*d/sqrt(3);
            ypart3 = 2*d/sqrt(3);
            ypart4 = 2*d/sqrt(3);
        else
            m = ny(j) - ny(i);
            ypart1 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(1/4-d/3)) + exp(1i*2*pi*(m)*(1/4+d/3)));
            ypart2 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(3/4-d/3)) + exp(1i*2*pi*(m)*(3/4+d/3)));
            ypart3 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(3/4-d/3)) + exp(1i*2*pi*(m)*(3/4+d/3)));
            ypart4 = 1./(2*pi*1i*(m)) * (-exp(1i*2*pi*(m)*(1/4-d/3)) + exp(1i*2*pi*(m)*(1/4+d/3)));
        end

        h(i,j) = h(i,j) + v * xpart1 * ypart1 + v * xpart2 * ypart2 + v * xpart3 * ypart3 + v * xpart4 * ypart4;
    end
end

% Fill the diagonals, taking the conjugate transpose for the lower diag
h = h + h';

% Rectangular K-space values
Kx = -sqrt(3)*pi:(sqrt(3)*pi/(iter)):sqrt(3)*pi;
Ky = -pi:(pi/(iter)):pi;

for i = 1:2*iter+1
    for j = 1:2*iter+1
        bloch = 4*(n(:,2).*n(:,2) + 3*n(:,3).*n(:,3)) + 16/(3*sqrt(3))*v*d*d + (4/pi)*(nx*Kx(i) + 3*ny*Ky(j)) + (Kx(i)*Kx(i) + 3*Ky(j)*Ky(j))/(pi*pi);
        Ebloch = eig (h + diag(bloch));
        E(i,j,:) = Ebloch;
    end
end

figure;
[X, Y] = meshgrid(Kx/pi, Ky/pi);
surf(Y,X,E(:,:,1)); hold on;
surf(Y,X,E(:,:,2));
surf(Y,X,E(:,:,3));
surf(Y,X,E(:,:,4));
zlabel('E/E_{ISW}', 'FontSize', 18);
ylabel('K_{y}a/\pi','FontSize', 16);
xlabel('K_{x}a/\pi','FontSize', 16);
end
```
