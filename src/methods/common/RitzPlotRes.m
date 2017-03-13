function RitzPlotRes(k, res, theta)
    % creates a socalled Ritz plot
    % k lists the order of the Krylov space
    % res is a cell array: res{i} contains the residuals of the Krylov 
    % space of order k(i)
    % theta is a cell array: theta{i} contains the Ritz values of the
    % Krylov space of order k(i)
    % daan.camps@cs.kuleuven.be
    t = [1e-1, 1e-3, 1e-6];
    bl = [0 0 153/256];
    gr = [0 102/256 0];
    or = [204/256 51/256 0];
    re = [204/256 0 0];
    if isreal(theta{1}),
        figure;
        for i=1:length(k)
            resc = res{i}; thetac = theta{i};
%             b = thetac(resc>t(1));
%             g = thetac(resc<=t(1)); %g = g(resc>t(2));
%             o = thetac(resc<=t(2)); %o = o(resc>t(3));
%             r = thetac(resc<=t(3));
            [b,g,o,r]= subdivideData(thetac,resc,t);
            plot(k(i)*ones(size(b)),b,'x','Color',bl);
            if (~isempty(b) && i==1), hold on; end
            plot(k(i)*ones(size(g)),g,'x','Color',gr);
            if (isempty(b) && ~isempty(g) && i==1), hold on; end
            plot(k(i)*ones(size(o)),o,'x','Color',or);
            if (isempty(b) && isempty(g) && ~isempty(o) && i==1), hold on; end
            plot(k(i)*ones(size(r)),r,'x','Color',re);
        end
        hold off;
    else
        figure;
        for i=1:length(k)
            resc = res{i}; thetac = theta{i};
%             b = thetac(resc>t(1));
%             g = thetac(resc<=t(1)); %g = g(resc>t(2));
%             o = thetac(resc<=t(2)); %o = o(resc>t(3));
%             r = thetac(resc<=t(3));
            [b,g,o,r]= subdivideData(thetac,resc,t);
            plot3(real(b),imag(b),k(i)*ones(size(b)),'x','Color',bl);
            if (~isempty(b) && i==1), hold on; end
            plot3(real(g),imag(g),k(i)*ones(size(g)),'x','Color',gr);
            if (isempty(b) && ~isempty(g) && i==1), hold on; end
            plot3(real(o),imag(o),k(i)*ones(size(o)),'x','Color',or);
            if (isempty(b) && isempty(g) && ~isempty(o) && i==1), hold on; end
            plot3(real(r),imag(r),k(i)*ones(size(r)),'x','Color',re);
        end
        hold off;
    end
end

function [b,g,o,r]= subdivideData(theta,res,t)
    b = []; g = []; o = []; r = [];
    for i=1:length(res)
       if res(i)>t(1),
           b(end+1) = theta(i);
       elseif (res(i)<=t(1)) && (res(i)>t(2)),
           g(end+1) = theta(i);
       elseif (res(i)<=t(2)) && (res(i)>t(3)),
           o(end+1) = theta(i);
       else
           r(end+1) = theta(i);
       end
    end
end
