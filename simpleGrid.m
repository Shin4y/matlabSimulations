
a = 0.05;     % infectivity
b = 0.01;     % recovery rate
tmax = 180;
clockmax = 250;
dt = tmax / clockmax;


S = [4000;4000;4000;4000;4000;4000]
I = [1500;1500;1500;1500;1500;1500]
R = [2000;2000;2000;2000;2000;2000] 
vaccRate = [0; 0.1; 0.2; 0.3; 0.4; 0.5]

NodeTable = table(S,I,R, vaccRate);
startNodes = [1 1 1 2 2 2 3 4 6];
endNodes =   [2 3 4 3 4 5 4 5 1];
weights =    [2 2 1 1 2 3 1 3 4];


moveTo1 = 0;
G = graph(startNodes, endNodes, weights, NodeTable)

disp(G.Nodes.R(1))
disp(G.Nodes.S(1))
disp(G.Edges)
numberOfNodes = 6;

tsave = zeros(numberOfNodes, clockmax);
Ssave = zeros(numberOfNodes, clockmax);
Isave = zeros(numberOfNodes, clockmax);
Rsave = zeros(numberOfNodes, clockmax);
rng('shuffle')

%vidObj = VideoWriter('infect.avi');
%open(vidObj);

for clock = 1:clockmax
    t = clock * dt;
    totalWeight = 0;
    
    for node = 1: numberOfNodes
        lowerCheck = 0;
        upperCheck = G.Edges.Weight(1);
        done = false;
        
       
        
        %disp(G.Nodes.I(node));
        for i = 1:10:G.Nodes.I(node)
            destinationNode = -1;
            if(randi(20) < 5)
                destinationNode = whichOne(G, node, numberOfNodes);
            end
            
            
            
            if(destinationNode ~= -1 && G.Nodes.I(node) > 30)
                G.Nodes.I(node) = G.Nodes.I(node) - 10;
                G.Nodes.I(destinationNode) = G.Nodes.I(destinationNode) + 10;
                
            end
        end
        
        
        for rode = 1: numberOfNodes
            I = G.Nodes.I(rode);
            S = G.Nodes.S(rode);
            R = G.Nodes.R(rode);
            %disp(I);
            %disp(R);
            %disp(S);
            N = I + S + R;
            
            
            StoI = dt * a * I/N * S;
            
            ItoR = dt * b * I;
            
            S = S - StoI;
            
            I = I + StoI - ItoR;
            R = R + ItoR + vacc(G, node);
            G.Nodes.I(rode) = I;
            G.Nodes.S(rode) = S;
            G.Nodes.R(rode) = R;
            tsave(rode, clock) = t;
            Ssave(rode, clock) = S;
            Isave(rode, clock) = I;
            Rsave(rode, clock) = R;
            
            %G.Nodes.NodeColors = G.Nodes.I;
            %G.Edges.LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
            %p = plot(G);
            %p.NodeLabel = G.Nodes.I;
            %p.NodeCData = G.Nodes.NodeColors;
            %p.LineWidth = G.Edges.LWidths;

            

            %colorbar
        
            %currFrame = getframe;
            %writeVideo(vidObj,currFrame);
            %delete(p) 
            
        end
        
    end
    
end

close(vidObj);
disp(G.Nodes.S);

f1 = figure;
 
subplot(2, 3, 1)
p1 = plot(tsave(1,:), Isave(1,:), 'r');
hold on;
p2 = plot(tsave(1,:), Ssave(1,:), 'g');
p3 = plot(tsave(1,:), Rsave(1,:), 'b');
hold off;
title("node1")
 
linehandles = [p1, p2, p3];
cols = cell2mat(get(linehandles, 'color'));
[~, uidx] = unique(cols, 'rows', 'stable');
legend(linehandles(uidx), {'I', 'S', 'R'})
 
subplot(2, 3, 2)
p4 = plot(tsave(2,:), Isave(2,:), 'r');
hold on;
p5 = plot(tsave(2,:), Ssave(2,:), 'g');
p6 = plot(tsave(2,:), Rsave(2,:), 'b');
hold off;
title("node2")
 
 
linehandles = [p4, p5, p6];
cols = cell2mat(get(linehandles, 'color'));
[~, uidx] = unique(cols, 'rows', 'stable');
legend(linehandles(uidx), {'I', 'S', 'R'})
 
subplot(2, 3, 3)
p7 = plot(tsave(3,:), Isave(3,:), 'r');
hold on;
p8 = plot(tsave(3,:), Ssave(3,:), 'g');
p9 = plot(tsave(3,:), Rsave(3,:), 'b');
hold off;
title("node3")
 
linehandles = [p7, p8, p9];
cols = cell2mat(get(linehandles, 'color'));
[~, uidx] = unique(cols, 'rows', 'stable');
legend(linehandles(uidx), {'I', 'S', 'R'})
 
subplot(2, 3, 4)
p10 = plot(tsave(4,:), Isave(4,:), 'r');
hold on;
p11 = plot(tsave(4,:), Ssave(4,:), 'g');
p12 = plot(tsave(4,:), Rsave(4,:), 'b');
hold off;
title("node4")
 
linehandles = [p10, p11, p12];
cols = cell2mat(get(linehandles, 'color'));
[~, uidx] = unique(cols, 'rows', 'stable');
legend(linehandles(uidx), {'I', 'S', 'R'})
 
subplot(2, 3, 5)
p13 = plot(tsave(5,:), Isave(5,:), 'r');
hold on;
p14 = plot(tsave(5,:), Ssave(5,:), 'g');
p15 = plot(tsave(5,:), Rsave(5,:), 'b');
hold off;
title("node5")
 
linehandles = [p13, p14, p15];
cols = cell2mat(get(linehandles, 'color'));
[~, uidx] = unique(cols, 'rows', 'stable');
legend(linehandles(uidx), {'I', 'S', 'R'})
 
subplot(2, 3, 6)
p16 = plot(tsave(6,:), Isave(6,:), 'r');
hold on;
p17 = plot(tsave(6,:), Ssave(6,:), 'g');
p18 = plot(tsave(6,:), Rsave(6,:), 'b');
hold off;
title("node6")
 
linehandles = [p16, p17, p18];
cols = cell2mat(get(linehandles, 'color'));
[~, uidx] = unique(cols, 'rows', 'stable');
legend(linehandles(uidx), {'I', 'S', 'R'})
 
savefig(f1, "meme");



function ret = whichOne(G, node, totalNodes)
rng('shuffle')
A = adjacency(G);
totalWeight = 0;
for i = 1:totalNodes
    
    if (A(node, i) == 1)
      
        totalWeight = totalWeight + G.Edges.Weight(findedge(G, node, i));
    end
end


index = 1;
storage = zeros(1, totalWeight);
for s = 1:totalNodes
    if (A(node, s) == 1)
        for r = 1:G.Edges.Weight(findedge(G, node, s)) + 1
            storage(index) = s;
            index = index + 1;
        end
    end
end

choice = randi(size(storage));
%disp(storage)
%formatSpec = "choice is: %d. ret is: %d. maxchoice is: %d. index is %d\n";
%fprintf(formatSpec, choice, storage(choice), totalWeight, index);

ret = storage(choice);
end


function z = vacc(G, node)
    z = (G.Nodes.S(node)/50) * G.Nodes.vaccRate(node);
end







