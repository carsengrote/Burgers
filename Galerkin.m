L = 2*pi; % Length of the interval, [0,L]
N = 100; % Number of points spatially
h = L/N; % Space discretization to evaluate function at
% Fundamental frequency is 2pi/L
x = 0:h:L - h; % Discretization of interval

createVideo = 1;
if (createVideo == 1)
    writerObj = VideoWriter('BurgersVideo.mp4','MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = 20;
    open(writerObj);
end

T = 1; % Final time
k = T/30000; % Time step, 4000 works with 100 points
% Initial conditions at t = 0
height = 2; % Initial wave height (height of 3, viscosity of 5 is a fun near blow up)
center = pi-1;
u = exp(-1*(center-x).^2)*height; % Gauss wave centered at x = center, with given height 
v = 2; % Viscosity
a = fft(u); % Get the fourier coeffficients, if we have u at N points
            % then we get the a_k for -N/2 <= k <= N/2 - 1

uPrime = zeros(1,length(a)); % Placeholder for the derivative values

for n = 0: T/k - 1

    uPrime = Derive(a,uPrime,v);
    a = a + k*uPrime; % Euler's Timestep

    if mod(n,4) == 0 % Can change the modulus number to make the video faster or slower

        if createVideo == 1
            tiledlayout(2,1);
            nexttile;
            plot(x, real(ifft(a)),'LineWidth',2)
            xlabel('$$x$$','interpreter','latex');
            ylabel('$$u(x,t)$$','interpreter','latex');
            xlim([0 (L-h)]);
            ylim([-1 height + .2]);
            set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right');
            title(['t = ', num2str(round((n*k)*100)/100) 's,  v = ',num2str(v)]);

            nexttile;
            plot(x(2:N/4+1)./h, abs(a(2:N/4+1)), 'bo','MarkerFaceColor','b')
            xlabel('$$k$$','interpreter','latex');
            ylabel('$$|a_{k}|$$','interpreter','latex');
            xlim([0,N/4 + 1])
            ylim([0 10]);
            set(get(gca,'ylabel'),'rotation',0,'HorizontalAlignment','right');
            set(findall(gcf,'-property','FontSize'),'FontSize',25);
            set(gcf,'position',[100,100,1000,1000])
            writeVideo(writerObj, getframe(gcf));
        end   
    end
end

if createVideo == 1
     close(writerObj);
end