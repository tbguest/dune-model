
% Werner1995.m


clear
close all

set(0, 'defaultaxesfontsize', 14)

% for iii = 1:3

% lattice
delx = 1;
xmax = 999;
x = 0:delx:xmax;
nx = length(x);

h = zeros(1, nx);

SAR = 1/3; % slab aspect ratio (< of repose = atan(2/3))
mean_slabs = 1; % mean number of slabs/lattice point

Ps = 0.6;
Pns = 0.4;
L = 5;

tmax = 400;

%figure(10), clf, hold on


% Setup step:
% populate lattice with randomly placed slabs
for ii = 1:mean_slabs*nx
    Ir = round(1 + (nx-1)*rand(1,1));
    h(Ir) = h(Ir) + 1;
    
    % check angle of repose
    stable = 0; % conditions not yet met
    while stable == 0
        if isequal(Ir, 1) 
            Ilo = nx;
            Ihi = 2;
        elseif isequal(Ir, nx)
            Ilo = nx-1;
            Ihi = 1;
        else
            Ilo = Ir - 1;
            Ihi = Ir + 1;
        end
        d1 = h(Ir) - h(Ilo);
        d2 = h(Ir) - h(Ihi);
        if abs(d1) > 1 && abs(d2) > 1
            h(Ir) = h(Ir) - 1;
            if isequal(d1, d2)
                Ipickone = round(1 + rand(1,1));
                if isequal(Ipickone, 1)
                    h(Ilo) = h(Ilo) + 1;
                    Ir = Ilo;
                else
                    h(Ihi) = h(Ihi) + 1;
                    Ir = Ihi;
                end
            end
        elseif abs(d1) > 1 || abs(d2) > 1
            h(Ir) = h(Ir) - 1;
            if abs(d1) > abs(d2)
                h(Ilo) = h(Ilo) + 1;
                Ir = Ilo;
            else
                h(Ihi) = h(Ihi) + 1;
                Ir = Ihi;
            end
        % condition has been met    
        else
            stable = 1;
        end
    end % while 1
         
end % setup step

% figure(1), clf
% plot(x, h)
% colorbar
% 
% figure(2), clf
% plot(diff(h))


%%%%%%%%% main model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tcount = 0;
allcount = 0;

while tcount < tmax
    
    % saltation step:     
    Ir = round(1 + (nx-1)*rand(1,1));
    
    % if the substrate is not exposed:
    if h(Ir) > 0
        
        % check if in shadow zone:
        Icheck = Ir - [1:nx-2];
        jnk = find(Icheck<1);
        for trsh = 1:length(jnk)
            Icheck(jnk(trsh)) = Icheck(jnk(trsh)) + nx-1;
        end
        
        counter = 0;
        for jj = 1:length(Icheck)
            if atan(h(Icheck(jj))/(jj*3/2)) > 15*pi/180
                counter = counter + 1;
            end
        end
        
        if counter == 0 % not in shadow zone; proceed
        
            h(Ir) = h(Ir) - 1;    

            % check angle of repose:
            stable = 0; % conditions not yet met
            while stable == 0
                if isequal(Ir, 1) 
                    Ilo = nx;
                    Ihi = 2;
                elseif isequal(Ir, nx)
                    Ilo = nx-1;
                    Ihi = 1;
                else
                    Ilo = Ir - 1;
                    Ihi = Ir + 1;
                end
                d1 = h(Ir) - h(Ilo);
                d2 = h(Ir) - h(Ihi);
                if abs(d1) > 1 && abs(d2) > 1
                    h(Ir) = h(Ir) + 1;
                    if isequal(d1, d2)
                        Ipickone = round(1 + rand(1,1));
                        if isequal(Ipickone, 1)
                            h(Ilo) = h(Ilo) - 1;
                            Ir = Ilo;
                        else
                            h(Ihi) = h(Ihi) - 1;
                            Ir = Ihi;
                        end
                    end
                elseif abs(d1) > 1 || abs(d2) > 1
                    h(Ir) = h(Ir) + 1;
                    if abs(d1) > abs(d2)
                        h(Ilo) = h(Ilo) - 1;
                        Ir = Ilo;
                    else
                        h(Ihi) = h(Ihi) - 1;
                        Ir = Ihi;
                    end
                % condition has been met    
                else
                    stable = 1;
                end
            end % while 1

            % grain saltates:
            transport = 1;

            while transport == 1
                Ir = Ir + L;
                if Ir > nx
                    Ir = mod(Ir, nx);
                end

                % prob of deposition
                
                %%%%%%%%%%%%%% first, check if in shadow zone%%%%%%
                % check if in shadow zone:
                Icheck = Ir - [1:nx-2];
                jnk = find(Icheck<1);
                for trsh = 1:length(jnk)
                    Icheck(jnk(trsh)) = Icheck(jnk(trsh)) + nx-1;
                end

                counter = 0;
                for jj = 1:length(Icheck)
                    if atan(h(Icheck(jj))/(jj*3/2)) > 15*180/pi
                        counter = counter + 1;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if counter > 0
                    h(Ir) = h(Ir) + 1; % deposited with P=1
                        transport = 0;
                
                elseif h(Ir) > 0
                    if rand(1,1) < Ps
                        h(Ir) = h(Ir) + 1;
                        transport = 0;
                    end
                else
                    if rand(1,1) < Pns
                        h(Ir) = h(Ir) + 1;
                        transport = 0;
                    end  
                end
            end % transport

            % check angle of repose
            stable = 0; % conditions not yet met
            while stable == 0
                if isequal(Ir, 1) 
                    Ilo = nx;
                    Ihi = 2;
                elseif isequal(Ir, nx)
                    Ilo = nx-1;
                    Ihi = 1;
                else
                    Ilo = Ir - 1;
                    Ihi = Ir + 1;
                end
                d1 = h(Ir) - h(Ilo);
                d2 = h(Ir) - h(Ihi);
                if abs(d1) > 1 && abs(d2) > 1
                    h(Ir) = h(Ir) - 1;
                    if isequal(d1, d2)
                        Ipickone = round(1 + rand(1,1));
                        if isequal(Ipickone, 1)
                            h(Ilo) = h(Ilo) + 1;
                            Ir = Ilo;
                        else
                            h(Ihi) = h(Ihi) + 1;
                            Ir = Ihi;
                        end
                    end
                elseif abs(d1) > 1 || abs(d2) > 1
                    h(Ir) = h(Ir) - 1;
                    if abs(d1) > abs(d2)
                        h(Ilo) = h(Ilo) + 1;
                        Ir = Ilo;
                    else
                        h(Ihi) = h(Ihi) + 1;
                        Ir = Ihi;
                    end
                % condition has been met    
                else
                    stable = 1;
                end
            end % while 1: < of repose
                
        end % shadow zone
        
    end % if h(Ir) > 0
         
    tcount = tcount + 1/nx;
    allcount = allcount + 1;
    
    if abs(max(diff(h(:)))) > 1
        disp('repose condition not met')
        break
    end
    
    ind = 0;
    if mod((allcount), 50000) == 1
        figure(10), hold on
%         subplot(3, 1, iii), hold on
        plot(x, h + tcount/20-1, 'k', 'linewidth', 1.5)
        axis tight
        box on
%         if iii ~= 3
%             set(gca, 'xticklabel', [])
%         end
    end
    

    
end % 

% end % iii

% figure(1), 
% subplot(3, 1, ii), hold on
% plot(x, h)

figure(10)
box on

