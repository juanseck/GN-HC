%  Global Neighborhood with Hill Climbing Algorithm (GN-HC)
%
%  Source codes demo version 2.0
%
%  Developed in MATLAB R2015a(7.08)
%
%  Author and programmer: 
%  Juan Carlos Seck Tuoh Mora
%  Nayeli Jazmin Escamilla Serna
%
% email:   jseck@uaeh.edu.mx
%          juanseck@gmail.com
%
%
%  Main paper:
%  A global neighborhood with hill climbing algorithm for fuzzy flexible job shop scheduling problem
%  Mathematics MDPI, DOI: http://
%_______________________________________________________________________________________________
% The initial parameters that you need are:
%__________________________________________
% numIndividuos = number of smart-cells
% numGeneraciones = maximum number of iterations
% numEstancamiento = maximum number of stagnation iterations
% probElitista = probability of elitist smart-cells
% tamVecindad = number of neighbors
% iteracionesTotalesEscalada = hill climbing iterations;
% iteracionesReinicioEscalada = hill climbing restart;
%______________________________________________________________________________________________

clear all
rng shuffle

%Parameters
numIndividuos = 80;
numGeneraciones = 500;
numEstancamiento = 250;
probElitista = 0.05;
tamVecindad = 5;
iteracionesTotalesEscalada = 150;
iteracionesReinicioEscalada = 15;

%Instance to be solved
problema='Instancia1_FFJSSP.mat';

%Read instance data
[numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles] = leerDatosProblema(problema);
%GN-HC algorithm
[mejorSO, mejorSM, mejorMakespan, PoblacionSO, PoblacionSM, PoblacionMakespan, convergencia] = GN_HC(numIndividuos, numGeneraciones, numEstancamiento, probElitista, tamVecindad, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles, iteracionesTotalesEscalada, iteracionesReinicioEscalada, 1);
%Display result
disp(['Problem: ' problema '; Makespan:(' num2str(mejorMakespan(1)) ',' num2str(mejorMakespan(2)) ',' num2str(mejorMakespan(3)) ')' ])
diagramaDeGanttMaquinas(mejorSO, mejorSM, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorInicioOperaciones, tablaTiempos)



%Read instance data
function [numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles] = leerDatosProblema(nombreArchivo)
%   Data format:
%   number of machines
%   number of jobs
%   operations per job
%   time table
archivo=load(nombreArchivo);
numeroTrabajos=archivo.numero_trabajos;
numeroMaquinas=archivo.numero_maquinas;
vectorNumOperaciones=archivo.operaciones_por_trabajo;
vectorOperaciones=[];
for nt=1:numeroTrabajos
    vectorInicioOperaciones(nt)=sum(vectorNumOperaciones(1:nt-1));
    operacionesTrabajo=ones(1,vectorNumOperaciones(nt))*nt;
    vectorOperaciones=[vectorOperaciones operacionesTrabajo];
end
tablaTiempos=archivo.tabla_tiempos;
numeroOperaciones=length(vectorOperaciones);
tablaMaquinasFactibles=[];
for oper=1:length(tablaTiempos)
    %Feasible machines
    indices_factibles = tablaTiempos(oper,1:3:numeroMaquinas*3) ~= 0;
    tablaMaquinasFactibles=[tablaMaquinasFactibles; indices_factibles];
end
end

%GN_HC algorithm
function [mejorSO, mejorSM, mejorMakespan, PoblacionSO, PoblacionSM, PoblacionMakespan, convergencia] = GN_HC(numIndividuos, numGeneraciones, numEstancamiento, probElitista, tamVecindad, numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones, vectorInicioOperaciones, vectorOperaciones, tablaTiempos, tablaMaquinasFactibles,iteracionesTotalesEscalada,iteracionesReinicioEscalada,bandImp)

%Smart-cell array
PoblacionSO=zeros(numIndividuos,numOperaciones);
PoblacionSM=zeros(numIndividuos,numOperaciones);
PoblacionMakespan=zeros(numIndividuos,3);
%Jobs table,
%Columns ordered by job and operation  (J_11, J_12, ... Jnm-1, Jnm)
%Row order:
%Assigned machine
%Processing position in the assigned machine
%Minimum processing time
%Most probable processing time
%Maximum processing time
%Operation position in SO
%Position of the previous operation in the table
PoblacionTablaTrabajos=zeros(10,numOperaciones,numIndividuos);

%Machine table,
%Columns ordered by machine and operation (M_11, M_12, ... Mmo-1, Jmo)
%Row order:
%Job assigned to the machine
%Operation of the job
%Minimum processing time
%Most probable processing time
%Maximum processing time
%Operation position in SO
PoblacionTablaMaquinas=zeros(9,numOperaciones,numIndividuos);
PoblacionVectorMaquinas=zeros(numIndividuos,numeroMaquinas);
PoblacionVectorOrdenMaq=zeros(numIndividuos,numOperaciones);
PoblacionPosMk=zeros(numIndividuos,1);
PoblacionPosTT=zeros(numIndividuos,numOperaciones);
PoblacionPosTM=zeros(numIndividuos,numOperaciones);

%Random initial population of smart-cells
for i=1:numIndividuos
    [PoblacionSO(i,:),PoblacionSM(i,:)] = generarIndividuoAleatorio(numeroMaquinas,numOperaciones,vectorOperaciones,tablaMaquinasFactibles);
end

%Evaluate population
[PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = calificarPoblacion(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numeroTrabajos, numeroMaquinas, numOperaciones,numIndividuos,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, 1);
%Select the best smart-cell
[mejorSO, mejorSM, mejorMakespan] = mejorIndividuo(PoblacionSO, PoblacionSM, PoblacionMakespan);

%Elitist smart-cells
numIndEl=round(numIndividuos*probElitista);
if mod(numIndividuos-numIndEl,2)==1
    numIndEl = numIndEl+1;
end

contIt=1;
contEst=1;
banderaCiclo=1;

%Convergence vector
convergencia=[];
convergencia(contIt,:)=mejorMakespan;
 
%Optimization loop
while(banderaCiclo)
    %Selection
    [PoblacionSO, PoblacionSM, PoblacionMakespan,PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = seleccion(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numIndividuos, numIndEl, 2);
    %GN-HC neigborhood
    [PoblacionSO, PoblacionSM, PoblacionMakespan,PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = aplicacionVecindadGN_HC(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, tamVecindad, numIndividuos, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, numIndEl);
    %Hill climbing
    [PoblacionSO, PoblacionSM, PoblacionMakespan,PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = busquedaPoblacionalEscaladaColina(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, numIndividuos, iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, numIndEl);
    %Select best solution in the iteration
    [nuevaMejorSO, nuevaMejorSM, nuevaMejorMakespan] = mejorIndividuo(PoblacionSO, PoblacionSM, PoblacionMakespan);
    %If improved, update best solution 
    [~,pos] = mayor_difuso(mejorMakespan,nuevaMejorMakespan);
    if pos==1
        suma=sum(mejorMakespan==nuevaMejorMakespan);
        if suma<3
            mejorSO = nuevaMejorSO;
            mejorSM = nuevaMejorSM;
            mejorMakespan = nuevaMejorMakespan;
            contEst=1;
        else
            %Increase stagnation
            contEst=contEst+1;
        end
    else
        contEst=contEst+1;
    end
    %Halt condition
    if ((contEst>=numEstancamiento) || (contIt>=numGeneraciones))
        banderaCiclo=0;
    end
    %Print best fuzzy makespan
    if mod(contIt,20)==0 && bandImp==1
        disp(['Iteracion: ' num2str(contIt) ' Makespan: ' num2str(mejorMakespan)]) 
    end
    contIt=contIt+1;
    %Keep convergence
    convergencia(contIt,:)=mejorMakespan;
end

end

%Random smart-cell
function [so,sm] = generarIndividuoAleatorio(numeroMaquinas,numOperaciones,vectorOperaciones,tablaMaquinasFactibles)
indices=randperm(numOperaciones);
so=vectorOperaciones(indices);
lista_maquinas = 1:numeroMaquinas;
sm=zeros(1,numOperaciones);
for i=1:numOperaciones
    maquinas=tablaMaquinasFactibles(i,:);
    factibles = lista_maquinas(logical(maquinas));
    numMaquinas=length(factibles);
    indMaquina=randi([1 numMaquinas]);
    maquina=factibles(indMaquina);
    sm(i)=maquina;
end

end

%Evaluate population
function [PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM] = calificarPoblacion(PobSO, PobSM, PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM, numeroTrabajos, numeroMaquinas, numOperaciones,numIndividuos,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, inicio)
for i=inicio:numIndividuos
    [PobMk(i,:), PobPMk(i), PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:),PobVOM(i,:), PobPTT(i,:), PobPTM(i,:)] = calcularMakespanActivoTablas(PobSO(i,:),PobSM(i,:),numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:), PobVOM(i,:), PobPTT(i,:), PobPTM(i,:));
end
end

%Calculate fuzzy makespan
function [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanActivoTablas(so,sm,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM)
makespan = zeros(1,3);
posmk = 0;
tiempoActualTrabajo=zeros(numeroTrabajos,3);
tiempoActualMaquina=zeros(numeroMaquinas,3);
operActualTrabajo=ones(1,numeroTrabajos);
posOperActualTrabajo=zeros(1,numeroTrabajos);
operActualMaquina=zeros(1,numeroMaquinas);
posOperActualMaquina=zeros(1,numeroMaquinas);
Ttoms=zeros(numeroMaquinas,3*(numeroTrabajos*2));
Ttopms=zeros(numeroMaquinas,3*(numeroTrabajos*2));
for i=1:numeroMaquinas
    suma=sum(sm==i);
    vectorMaquinas(i)=suma;
end
for i=1:numOperaciones
    trabajo=so(i);
    operacion=operActualTrabajo(trabajo);
    operActualTrabajo(trabajo)=operActualTrabajo(trabajo)+1;
    posicion=vectorInicioOperaciones(trabajo)+operacion;
    maquina=sm(posicion);
    operActualMaquina(maquina)=operActualMaquina(maquina)+1;
    operacionMaquina=operActualMaquina(maquina);
    tiempo=tablaTiempos(posicion,(maquina-1)*3+1:(maquina-1)*3+3);
    aux1=tiempoActualTrabajo(trabajo,:);
    pos_aux1=posOperActualTrabajo(trabajo);
    toat=aux1;
    noms=operActualMaquina(maquina)-1;
    aux2=tiempoActualMaquina(maquina,:);
    pos_aux2=posOperActualMaquina(maquina);
    [tiempoInicial,flag]=mayor_difuso(aux1,aux2);
    tiempoFinal=tiempoInicial+tiempo;
    tiempoActualMaquina(maquina,:)=tiempoFinal;
    tiempoActualTrabajo(trabajo,:)=tiempoFinal;
    if flag==1
        pos_oper=pos_aux1;
    else
        pos_oper=pos_aux2;
    end
    posOperActualTrabajo(trabajo)=posicion;
    posOperActualMaquina(maquina)=posicion;
    [makespan, flag] = mayor_difuso(makespan,tiempoFinal);
    if flag==2
        posmk = posicion;
    end
    Ttoms(maquina,noms*3+1:noms*3+3)=tiempoFinal;
    Ttopms(maquina,noms*3+1:noms*3+3)=tiempoInicial;
    TablaTrabajos(1,posicion) = maquina;
    TablaTrabajos(2,posicion) = operacionMaquina;
    TablaTrabajos(3,posicion) = tiempoFinal(1);
    TablaTrabajos(4,posicion) = tiempoFinal(2);
    TablaTrabajos(5,posicion) = tiempoFinal(3);
    TablaTrabajos(6,posicion) = i;
    TablaTrabajos(7,posicion) = pos_oper;
    TablaTrabajos(9,posicion) = tiempo(3); %average duration of the operation
    TablaTrabajos(10,posicion) = tiempoInicial(3); %initial average time of the operation
    posicionMaquina=sum(vectorMaquinas(1:maquina-1))+operacionMaquina;
    TablaMaquinas(1,posicionMaquina) = trabajo;
    TablaMaquinas(2,posicionMaquina) = operacion;
    TablaMaquinas(3,posicionMaquina) = tiempoFinal(1);
    TablaMaquinas(4,posicionMaquina) = tiempoFinal(2);
    TablaMaquinas(5,posicionMaquina) = tiempoFinal(3);
    TablaMaquinas(6,posicionMaquina) = i;
    TablaMaquinas(8,posicionMaquina) = tiempo(3); %average duration of the operation
    TablaMaquinas(9,posicion) = tiempoInicial(3); %initial average time of the operation
    vectorOrdenMaq(i) = maquina;
    posTT(i)=posicion;
    posTM(i)=posicionMaquina;
end
%Tail time
for i=numOperaciones:-1:1
    postt=posTT(i);
    postm=posTM(i);
    maq=TablaTrabajos(1,postt);
    posm=TablaTrabajos(2,postt);
    trab=TablaMaquinas(1,postm);
    oper=TablaMaquinas(2,postm);
    if oper == vectorNumOperaciones(trab)
        tiempo_t=0;
    else
        tiempo_t=TablaTrabajos(8,postt+1)+TablaTrabajos(9,postt+1);
    end
    if posm==vectorMaquinas(maq)
        tiempo_m=0;
    else
        tiempo_m=TablaMaquinas(7,postm+1)+TablaMaquinas(8,postm+1);
    end
    tiempo_cola=max([tiempo_t tiempo_m]);
    TablaTrabajos(8,postt)=tiempo_cola;
    TablaMaquinas(7,postm)=tiempo_cola;
end
end

%Select best smart-cell
function [mejorSO, mejorSM, mejorMk, indice] = mejorIndividuo(PobSO, PobSM, PobMk)
num_ind=size(PobMk,1);
indice=1;
mejorMk=PobMk(1,:);
mejorSO=PobSO(1,:);
mejorSM=PobSM(1,:);
for i=2:num_ind
    z2=PobMk(i,:);
    [~,pos] = mayor_difuso(mejorMk,z2);
    if pos==1
        mejorMk = PobMk(i,:);
        mejorSO = PobSO(i,:);
        mejorSM = PobSM(i,:);
        indice = i;
    end
end
end

%Selection by elitism and tournement
function [nuevaPobSO, nuevaPobSM, nuevaPobMk, nuevaPobPosMk, nuevaPobTablaTrabajos, nuevaPobTablaMaquinas, nuevaPobVectorMaquinas, nuevaPobVectorOrdenMaq, nuevaPobPosTT, nuevaPobPosTM] = seleccion(PobSO, PobSM, PobMk, PobPosMk, PobTablaTrabajos, PobTablaMaquinas, PobVectorMaquinas, PobVectorOrdenMaq, PobPosTT, PobPosTM, numsoluciones, numElitista, numCompetidores)
[nuevaPobSO, nuevaPobSM, nuevaPobMk, nuevaPobPosMk, nuevaPobTablaTrabajos, nuevaPobTablaMaquinas, nuevaPobVectorMaquinas, nuevaPobVectorOrdenMaq, nuevaPobPosTT, nuevaPobPosTM] = seleccionElitista(PobSO, PobSM, PobMk, PobPosMk, PobTablaTrabajos, PobTablaMaquinas, PobVectorMaquinas, PobVectorOrdenMaq, PobPosTT, PobPosTM, numElitista);
ind=numElitista+1;
while (ind<=numsoluciones)
    indGanador = seleccionTorneo(PobMk, numsoluciones, numCompetidores);
    nuevaPobSO(ind,:)=PobSO(indGanador,:);
    nuevaPobSM(ind,:)=PobSM(indGanador,:);
    nuevaPobMk(ind)=PobMk(indGanador);
    nuevaPobPosMk(ind)=PobPosMk(indGanador);
    nuevaPobTablaTrabajos(:,:,ind)=PobTablaTrabajos(:,:,indGanador);
    nuevaPobTablaMaquinas(:,:,ind)=PobTablaMaquinas(:,:,indGanador);
    nuevaPobVectorMaquinas(ind,:)=PobVectorMaquinas(indGanador,:);
    nuevaPobVectorOrdenMaq(ind,:)=PobVectorOrdenMaq(indGanador,:);
    nuevaPobPosTT(ind,:)=PobPosTT(indGanador,:);
    nuevaPobPosTM(ind,:)=PobPosTM(indGanador,:);
    ind=ind+1;
end

end

%Elitist selection
function [nuevaPobSO, nuevaPobSM, nuevaPobMk, nuevaPobPosMk, nuevaPobTablaTrabajos, nuevaPobTablaMaquinas, nuevaPobVectorMaquinas, nuevaPobVectorOrdenMaq, nuevaPobPosTT, nuevaPobPosTM] = seleccionElitista(PobSO, PobSM, PobMk, PobPosMk, PobTablaTrabajos, PobTablaMaquinas, PobVectorMaquinas, PobVectorOrdenMaq, PobPosTT, PobPosTM, numElitista)
nuevaPobSO=PobSO;
nuevaPobSM=PobSM;
nuevaPobMk=PobMk;
nuevaPobPosMk=PobPosMk;
nuevaPobTablaTrabajos=PobTablaTrabajos;
nuevaPobTablaMaquinas=PobTablaMaquinas;
nuevaPobVectorMaquinas=PobVectorMaquinas;
nuevaPobVectorOrdenMaq=PobVectorOrdenMaq; 
nuevaPobPosTT=PobPosTT;
nuevaPobPosTM=PobPosTM;
[ indices ] = ordenarValoresDifusos(PobMk);
for i=1:numElitista
    nuevaPobSO(i,:)=PobSO(indices(i),:);
    nuevaPobSM(i,:)=PobSM(indices(i),:);
    nuevaPobMk(i)=PobMk(indices(i));
    nuevaPobPosMk(i)=PobPosMk(indices(i));
    nuevaPobTablaTrabajos(:,:,i)=PobTablaTrabajos(:,:,indices(i));
    nuevaPobTablaMaquinas(:,:,i)=PobTablaMaquinas(:,:,indices(i));
    nuevaPobVectorMaquinas(i,:)=PobVectorMaquinas(indices(i),:);
    nuevaPobVectorOrdenMaq(i,:)=PobVectorOrdenMaq(indices(i),:);
    nuevaPobPosTT(i,:)=PobPosTT(indices(i),:);
    nuevaPobPosTM(i,:)=PobPosTM(indices(i),:);
end
end

%Tournement selection
function indGanador= seleccionTorneo(PobMk,numsoluciones,numCompetidores)
indices=randperm(numsoluciones,numCompetidores);
z1=PobMk(indices(1),:);
z2=PobMk(indices(2),:);
[~,pos] = mayor_difuso(z1,z2);
ind=1;
if pos==1
    ind=2;
end
indGanador=indices(ind);
end

%Calculate GN neighborhood to all the population of smart-cells
function [PoblacionSO, PoblacionSM, PoblacionMakespan,PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM] = aplicacionVecindadGN_HC(PoblacionSO, PoblacionSM, PoblacionMakespan, PoblacionPosMk, PoblacionTablaTrabajos, PoblacionTablaMaquinas, PoblacionVectorMaquinas, PoblacionVectorOrdenMaq, PoblacionPosTT, PoblacionPosTM, tamVecindad, numIndividuos, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, inicio)
for i=inicio+1:numIndividuos
    so=PoblacionSO(i,:);
    sm=PoblacionSM(i,:);
    TablaTrabajos=PoblacionTablaTrabajos(:,:,i);
    TablaMaquinas=PoblacionTablaMaquinas(:,:,i);
    vectorMaquinas=PoblacionVectorMaquinas(i,:);
    vectorOrdenMaq=PoblacionVectorOrdenMaq(i,:);
    posTT=PoblacionPosTT(i,:);
    posTM=PoblacionPosTM(i,:);
    [mejor_so, mejor_sm, mejor_mk, mejor_posmk, mejor_TablaTrabajos, mejor_TablaMaquinas, mejor_vectorMaquinas, mejor_vectorOrdenMaq, mejor_posTT, mejor_posTM] = vecindadGN(so, sm, tamVecindad, PoblacionSO, numIndividuos, numeroTrabajos, numeroMaquinas, numeroOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
    PoblacionSO(i,:)=mejor_so;
    PoblacionSM(i,:)=mejor_sm;
    PoblacionMakespan(i,:)=mejor_mk;
    PoblacionPosMk(i)=mejor_posmk; 
    PoblacionTablaTrabajos(:,:,i)=mejor_TablaTrabajos; 
    PoblacionTablaMaquinas(:,:,i)=mejor_TablaMaquinas; 
    PoblacionVectorMaquinas(i,:)=mejor_vectorMaquinas; 
    PoblacionVectorOrdenMaq(i,:)=mejor_vectorOrdenMaq; 
    PoblacionPosTT(i,:)=mejor_posTT; 
    PoblacionPosTM(i,:)=mejor_posTM;
end
end

%Calculate GN neighborhood
function [mejor_so, mejor_sm, mejor_mk, mejor_posmk, mejor_TablaTrabajos, mejor_TablaMaquinas, mejor_vectorMaquinas, mejor_vectorOrdenMaq, mejor_posTT, mejor_posTM ] = vecindadGN(so, sm, tamVecindad, PoblacionSO, numIndividuos, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM)
mejor_mk=zeros(1,3);
for i=1:tamVecindad    
    opcion=rand;
    if opcion <= 0.5
        [nuevo_so, nuevo_sm, nuevo_mk, nuevo_posmk, nuevo_TablaTrabajos, nuevo_TablaMaquinas, nuevo_vectorMaquinas, nuevo_vectorOrdenMaq, nuevo_posTT, nuevo_posTM ] = insercionClasica(so, sm, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles,vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
    elseif opcion > 0.5 && opcion <= 0.75  
        [nuevo_so, nuevo_sm, nuevo_mk, nuevo_posmk, nuevo_TablaTrabajos, nuevo_TablaMaquinas, nuevo_vectorMaquinas, nuevo_vectorOrdenMaq, nuevo_posTT, nuevo_posTM ] = intercambioClasico(so, sm, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles,vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
    else
        ind=randi(numIndividuos);
        sof=PoblacionSO(ind,:);
        [nuevo_so, nuevo_sm, nuevo_mk, nuevo_posmk, nuevo_TablaTrabajos, nuevo_TablaMaquinas, nuevo_vectorMaquinas, nuevo_vectorOrdenMaq, nuevo_posTT, nuevo_posTM] = pathRelinkingClasico(so, sm, sof, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos,tablaMaquinasFactibles,vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
    end
    flag=0;
    if i==1
        flag=1;
    else
        [~,pos]=mayor_difuso(mejor_mk,nuevo_mk);
        if pos==1
            flag=1;
        end
    end
    if flag == 1
        mejor_so=nuevo_so;
        mejor_sm=nuevo_sm;
        mejor_mk=nuevo_mk;
        mejor_posmk=nuevo_posmk;
        mejor_TablaTrabajos=nuevo_TablaTrabajos; 
        mejor_TablaMaquinas=nuevo_TablaMaquinas; 
        mejor_vectorMaquinas=nuevo_vectorMaquinas; 
        mejor_vectorOrdenMaq=nuevo_vectorOrdenMaq; 
        mejor_posTT=nuevo_posTT; 
        mejor_posTM=nuevo_posTM;
    end
end

end

%Path-relinking
function [nuevo_so, nuevo_sm, nuevo_mk, nuevo_posmk, nuevo_TablaTrabajos, nuevo_TablaMaquinas, nuevo_vectorMaquinas, nuevo_vectorOrdenMaq, nuevo_posTT, nuevo_posTM] = pathRelinkingClasico(so, sm, sof, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos,tablaMaquinasFactibles,vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM)
posicionPR=randi(numOperaciones-0);
secuencia1=sof;
secuencia2=so;
secuenciaf = pathRelinkingSO(secuencia1,secuencia2,numOperaciones,posicionPR);
nuevo_so=secuenciaf;
nuevo_sm=sm;
if rand<=0.1
    [nuevo_sm] = mutacionMaquinas(numeroMaquinas, numOperaciones,sm,tablaMaquinasFactibles);
end
[nuevo_mk,nuevo_posmk,nuevo_TablaTrabajos,nuevo_TablaMaquinas,nuevo_vectorMaquinas,nuevo_vectorOrdenMaq,nuevo_posTT,nuevo_posTM] = calcularMakespanActivoTablas(nuevo_so,nuevo_sm,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
end

%Path -relinking over SO
function secuenciaf = pathRelinkingSO(secuencia1,secuencia2,longitudSec,posicion)
secuenciaf=secuencia1;
for i=1:posicion
    operacion=secuencia2(i);
    if secuenciaf(i)~=operacion
        j=i+1;
        if j>=longitudSec
            break
        end
        while(1)
            if secuenciaf(j)==operacion
                secuenciaf(j)=secuenciaf(i);
                secuenciaf(i)=operacion;
                break
            end
            j=j+1;
            if j>=longitudSec
                break
            end
        end
    end
end

end

%Insertion
function [nuevo_so, nuevo_sm, nuevo_mk, nuevo_posmk, nuevo_TablaTrabajos, nuevo_TablaMaquinas, nuevo_vectorMaquinas, nuevo_vectorOrdenMaq, nuevo_posTT, nuevo_posTM] = insercionClasica(so, sm, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles,vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM)
nuevo_so=so;
nuevo_sm=sm;
operacion=randi(numOperaciones);
ind=randi(numOperaciones);
while(ind==operacion)
    ind=randi(numOperaciones);
end
nuevo_so(ind)=so(operacion);
if operacion>ind
    nuevo_so(ind+1:operacion)=so(ind:operacion-1);
else
    nuevo_so(operacion:ind-1)=so(operacion+1:ind);
end
if rand<=0.1
    [nuevo_sm] = mutacionMaquinas(numeroMaquinas, numOperaciones,sm,tablaMaquinasFactibles);
end
[nuevo_mk,nuevo_posmk,nuevo_TablaTrabajos,nuevo_TablaMaquinas,nuevo_vectorMaquinas,nuevo_vectorOrdenMaq,nuevo_posTT,nuevo_posTM] = calcularMakespanActivoTablas(nuevo_so,nuevo_sm,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
end

%Swapping
function [nuevo_so, nuevo_sm, nuevo_mk, nuevo_posmk, nuevo_TablaTrabajos, nuevo_TablaMaquinas, nuevo_vectorMaquinas, nuevo_vectorOrdenMaq, nuevo_posTT, nuevo_posTM] = intercambioClasico(so, sm, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles,vectorNumOperaciones,TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM)
nuevo_so=so;
nuevo_sm=sm;
for i=1:8
    ind=randi(numOperaciones);
    operacion=nuevo_so(ind);
    ind2=randi(numOperaciones);
    operacion2=nuevo_so(ind2);
    while(operacion2==operacion)
        ind2=randi(numOperaciones);
        operacion2=nuevo_so(ind2);
    end
    nuevo_so(ind)=operacion2;
    nuevo_so(ind2)=operacion;
end
if rand<=0.1
    [nuevo_sm] = mutacionMaquinas(numeroMaquinas, numOperaciones,sm,tablaMaquinasFactibles);
end
[nuevo_mk,nuevo_posmk,nuevo_TablaTrabajos,nuevo_TablaMaquinas,nuevo_vectorMaquinas,nuevo_vectorOrdenMaq,nuevo_posTT,nuevo_posTM] = calcularMakespanActivoTablas(nuevo_so,nuevo_sm,numeroTrabajos, numeroMaquinas, numOperaciones, vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, TablaTrabajos, TablaMaquinas, vectorMaquinas, vectorOrdenMaq, posTT, posTM);
end

%Machine mutation
function [nuevaSol] = mutacionMaquinas(numeroMaquinas, numOper,solucion,tablaMaquinasFactibles)
nuevaSol = solucion;
lista_maquinas = 1:numeroMaquinas;
indices=randperm(numOper,floor(numOper/2));
for i=1:length(indices)
    maquinas=tablaMaquinasFactibles(indices(i),:);
    factibles = lista_maquinas(logical(maquinas));
    numMaquinas=length(factibles);
    indMaquina=randi([1 numMaquinas]);
    maq=factibles(indMaquina);
    nuevaSol(indices(i))=maq;
end

end

%Hill climbing for all the population of smart-cells
function [PobSO, PobSM, PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM] = busquedaPoblacionalEscaladaColina(PobSO, PobSM, PobMk, PobPMk, PobTT, PobTM, PobVM, PobVOM, PobPTT, PobPTM, numIndividuos, iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones, inicio)
for i=inicio+1:numIndividuos
    [PobSO(i,:), PobSM(i,:), PobMk(i,:), PobPMk(i), PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:), PobVOM(i,:), PobPTT(i,:), PobPTM(i,:)] = busquedaEscaladaColinaReinicioAleatorio(PobSO(i,:), PobSM(i,:), PobMk(i,:), PobPMk(i), PobTT(:,:,i), PobTM(:,:,i), PobVM(i,:), PobVOM(i,:), PobPTT(i,:), PobPTM(i,:), iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones);
end

end

%Hill climbing
function [mejorsolucionSO, mejorsolucionSM, mejorsolucionMk, mejorsolucionPMk, mejorsolucionTT, mejorsolucionTM, mejorsolucionVM, mejorsolucionVOM, mejorsolucionPTT, mejorsolucionPTM] = busquedaEscaladaColinaReinicioAleatorio(solucionSO, solucionSM, solucionMk, solucionPMk, solucionTT, solucionTM, solucionVM, solucionVOM, solucionPTT, solucionPTM, iteracionesTotalesEscalada, iteracionesReinicioEscalada, numeroTrabajos, numeroMaquinas, numOperaciones,vectorInicioOperaciones, tablaTiempos, tablaMaquinasFactibles, vectorNumOperaciones)
mejorsolucionSO = solucionSO;
mejorsolucionSM = solucionSM;
mejorsolucionMk = solucionMk;
mejorsolucionPMk = solucionPMk;
mejorsolucionTT = solucionTT;
mejorsolucionTM = solucionTM;
mejorsolucionVM = solucionVM;
mejorsolucionVOM = solucionVOM;
mejorsolucionPTT = solucionPTT;
mejorsolucionPTM = solucionPTM;
%Critical path
[RC, maqRC, tamRC] = rutaCriticaConTablas(numOperaciones, solucionTT, solucionPMk);
indSMRC = RC(tamRC:-1:1);
nuevas_SM=zeros(iteracionesReinicioEscalada,numOperaciones);
%Hill climbing loop
reinicio=0;
for i=1:iteracionesTotalesEscalada
    pos_ruleta=randi([1,tamRC]);
    pos_nueva=RC(pos_ruleta);
    maquina_actual=maqRC(pos_ruleta);
    nueva_maquina=maqRC(pos_ruleta);
    while(nueva_maquina==maquina_actual)
        indice=randi([1, numeroMaquinas]);
        if tablaMaquinasFactibles(pos_nueva,indice)==1
            nueva_maquina=indice;
        end
    end
    smAux=solucionSM;
    smAux(pos_nueva)=nueva_maquina;
    soAux=mejorsolucionSO;
    pos_os = solucionTT(6,pos_nueva);
    trabajo = soAux(pos_os);
    %Makespan estimation
    [nuevo_makespan] = estimarMakespan(pos_os,trabajo,nueva_maquina,vectorNumOperaciones,tablaTiempos,solucionTT,solucionTM,solucionPTT,solucionPTM,solucionVM,solucionVOM);
    if nuevo_makespan < mejorsolucionMk(3)        
        [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanActivoTablas(mejorsolucionSO,smAux,numeroTrabajos,numeroMaquinas,numOperaciones,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, solucionTT, solucionTM, solucionVM, solucionVOM, solucionPTT, solucionPTM);
        
        [~,pos] = mayor_difuso(makespan,mejorsolucionMk);
        if pos == 2
            mejorsolucionSM = smAux;
            mejorsolucionMk = makespan;
            mejorsolucionPMk = posmk;
            mejorsolucionTT = TablaTrabajos;
            mejorsolucionTM = TablaMaquinas;
            mejorsolucionVM = vectorMaquinas;
            mejorsolucionVOM = vectorOrdenMaq;
            mejorsolucionPTT = posTT;
            mejorsolucionPTM = posTM;
            solucionSM = smAux;
            solucionMk = makespan;
            solucionPMk = posmk;
            solucionTT = TablaTrabajos;
            solucionTM = TablaMaquinas;
            solucionVM = vectorMaquinas;
            solucionVOM = vectorOrdenMaq;
            solucionPTT = posTT;
            solucionPTM = posTM;
            reinicio=0;
            [RC, maqRC, tamRC] = rutaCriticaConTablas(numOperaciones, solucionTT, solucionPMk);
            indSMRC = RC(tamRC:-1:1);
        else
            reinicio = reinicio + 1;
            nuevas_SM(reinicio,:) = smAux;
        end
        %Restart hill climbing
        if reinicio >= iteracionesReinicioEscalada
            pos_al3=randi([1 reinicio]);
            smAux = nuevas_SM(pos_al3,:);
            [makespan,posmk,TablaTrabajos,TablaMaquinas,vectorMaquinas,vectorOrdenMaq,posTT,posTM] = calcularMakespanActivoTablas(mejorsolucionSO,smAux,numeroTrabajos,numeroMaquinas,numOperaciones,vectorNumOperaciones,vectorInicioOperaciones, tablaTiempos, mejorsolucionTT, mejorsolucionTM, mejorsolucionVM, mejorsolucionVOM, mejorsolucionPTT, mejorsolucionPTM);
            solucionSM = smAux;
            solucionMk = makespan;
            solucionPMk = posmk;
            solucionTT = TablaTrabajos;
            solucionTM = TablaMaquinas;
            solucionVM = vectorMaquinas;
            solucionVOM = vectorOrdenMaq;
            solucionPTT = posTT;
            solucionPTM = posTM;
            [RC, maqRC, tamRC] = rutaCriticaConTablas(numOperaciones, solucionTT, solucionPMk);
            indSMRC = RC(tamRC:-1:1);
            reinicio = 0;
        end
    end
end
end

%Critical path
function [RC, maqRC, tamRC] = rutaCriticaConTablas(numOperaciones, TablaTT, PosMk)
RC=zeros(1,numOperaciones);
maqRC=zeros(1,numOperaciones);
tamRC=0;
pos_crit = PosMk;
while(pos_crit > 0)
    tamRC=tamRC+1;
    RC(tamRC)=pos_crit;
    maqRC(tamRC)=TablaTT(1,pos_crit);
    pos_crit=TablaTT(7,pos_crit);
end
end

%Estimate maximum makespan
function [nuevo_makespan] = estimarMakespan(pos_os,trabajo,nueva_maq,vectorNumOperaciones,tablaTiempos,solucionTT,solucionTM,solucionPTT,solucionPTM,solucionVM,solucionVOM)
    ind_m=sum(solucionVOM(1:pos_os-1)==nueva_maq);
    ind_TM = sum(solucionVM(1:nueva_maq-1));
    if ind_m == 0
        tf_ma = 0;
    else
        ind_TM = ind_TM + ind_m;
        tf_ma = solucionTM(5,ind_TM);
    end
    ind_o = solucionTM(2,solucionPTM(pos_os));
    ind_TT = solucionPTT(pos_os);
    if ind_o == 1
        tf_oa = 0;
    else
        tf_oa = solucionTT(5,ind_TT-1);
    end
    max_tf_a = max([tf_ma, tf_oa]);    
    duracion = tablaTiempos(ind_TT,nueva_maq);
    if ind_m == solucionVM(nueva_maq)
        tc_ms = 0;
    else
        duracion_ms = solucionTM(8,ind_TM+1);
        tc_ms = duracion_ms + solucionTM(7,ind_TM+1);
    end
    if ind_o == vectorNumOperaciones(trabajo)
        tc_os = 0;
    else
        duracion_os = solucionTT(9,ind_TT+1);
        tc_os = duracion_os + solucionTT(8,ind_TT+1);
    end
    max_tf_s = max([tc_ms, tc_os]);
    nuevo_makespan = max_tf_a + duracion + max_tf_s;
end


%Fuzzy Gantt chart
function diagramaDeGanttMaquinas(so, sm, numeroTrabajos, numeroMaquinas, numOperaciones, vectorInicioOperaciones, tablaTiempos)
tiemposMinMaqAnt=zeros(numeroMaquinas,numeroTrabajos*2);
tiemposPromMaqAnt=zeros(numeroMaquinas,numeroTrabajos*2);
tiemposMaxMaqAnt=zeros(numeroMaquinas,numeroTrabajos*2);
tiemposMinMaq=zeros(numeroMaquinas,numeroTrabajos*2);
tiemposPromMaq=zeros(numeroMaquinas,numeroTrabajos*2);
tiemposMaxMaq=zeros(numeroMaquinas,numeroTrabajos*2);
trabMaq=zeros(numeroMaquinas,numeroTrabajos*2);
operTrabMaq=zeros(numeroMaquinas,numeroTrabajos*2);
makespan = zeros(1,3);
tiempoActualMaquina=zeros(3,numeroMaquinas);
operActualMaquina=zeros(1,numeroMaquinas);
tiempoActualTrabajo=zeros(3,numeroTrabajos);
operActualTrabajo=zeros(1,numeroTrabajos);
for i=1:numOperaciones
    trabajo=so(i);
    operActualTrabajo(trabajo)=operActualTrabajo(trabajo)+1;
    operacion=operActualTrabajo(trabajo);
    posicion=vectorInicioOperaciones(trabajo)+operacion;
    maquina=sm(posicion);
    tiempo=tablaTiempos(posicion,((maquina-1)*3)+1:((maquina-1)*3)+3);
    tiempoInicial=tiempoActualTrabajo(:,trabajo)';
    aux=tiempoActualMaquina(:,maquina)';
    [~,pos] = mayor_difuso(tiempoInicial,aux);
    if pos == 2
        tiempoInicial=aux;
    end
    tiempoFinal=tiempoInicial+tiempo;
    tiempoActualMaquina(:,maquina)=tiempoFinal';
    tiempoActualTrabajo(:,trabajo)=tiempoFinal';
    operActualMaquina(maquina)=operActualMaquina(maquina)+1;
    [~,pos] = mayor_difuso(makespan,tiempoFinal);
    if pos == 2
        makespan = tiempoFinal;
    end
    tiemposMinMaqAnt(maquina,operActualMaquina(maquina))=tiempoInicial(1);
    tiemposPromMaqAnt(maquina,operActualMaquina(maquina))=tiempoInicial(2);
    tiemposMaxMaqAnt(maquina,operActualMaquina(maquina))=tiempoInicial(3);
    tiemposMinMaq(maquina,operActualMaquina(maquina))=tiempoFinal(1);
    tiemposPromMaq(maquina,operActualMaquina(maquina))=tiempoFinal(2);
    tiemposMaxMaq(maquina,operActualMaquina(maquina))=tiempoFinal(3);
    trabMaq(maquina,operActualMaquina(maquina))=trabajo;
    operTrabMaq(maquina,operActualMaquina(maquina))=operActualTrabajo(trabajo);
end
%Graph
colorT=colormap(hsv(numeroTrabajos));
colorM=colormap(lines(numeroMaquinas));
figure(1);
clf
h=2.5;
for maquina=1:numeroMaquinas
    y=[maquina*h-0.75 maquina*h-0.25 maquina*h-0.25 maquina*h-0.75];
    x=[-4 -4 -1 -1 ];
    patch(x,y,colorM(maquina,:));
    etiqueta=strcat( 'M', num2str(maquina));
    text(-4.5,maquina*h-0.5,etiqueta)
    x=[0 makespan(3)+5];
    y=[maquina*h-0.5 maquina*h-0.5];
    line(x,y,'Color','k')
    for oper=1:operActualMaquina(maquina)
        x = [tiemposMinMaq(maquina,oper) tiemposPromMaq(maquina,oper) tiemposMaxMaq(maquina,oper)];
        y=[maquina*h-0.5 maquina*h-1.25 maquina*h-0.5];
        patch(x,y,colorT(trabMaq(maquina,oper),:));
        etiqueta=strcat( 'O_{', num2str(trabMaq(maquina,oper)), ',', num2str(operTrabMaq(maquina,oper)),'}' ); 
        text(tiemposPromMaq(maquina,oper)-2.5,maquina*h-1.5,etiqueta)
        etiqueta=strcat( '(', num2str(tiemposMinMaq(maquina,oper)), ',', num2str(tiemposPromMaq(maquina,oper)), ',', num2str(tiemposMaxMaq(maquina,oper)),')' );
        text(tiemposPromMaq(maquina,oper)-3,maquina*h-0.75,etiqueta)
        if tiemposMaxMaqAnt(maquina,oper)>0
            x = [tiemposMinMaqAnt(maquina,oper) tiemposPromMaqAnt(maquina,oper) tiemposMaxMaqAnt(maquina,oper)];
            y=[maquina*h-0.5 maquina*h+0.25 maquina*h-0.5];
            patch(x,y,colorT(trabMaq(maquina,oper),:));
            etiqueta=strcat( 'O_{', num2str(trabMaq(maquina,oper)), ',', num2str(operTrabMaq(maquina,oper)),'}' );
            text(tiemposPromMaqAnt(maquina,oper)-2.5,maquina*h+0.45,etiqueta)
            etiqueta=strcat( '(', num2str(tiemposMinMaqAnt(maquina,oper)), ',', num2str(tiemposPromMaqAnt(maquina,oper)), ',', num2str(tiemposMaxMaqAnt(maquina,oper)),')' );
            text(tiemposPromMaqAnt(maquina,oper)-3,maquina*h-0.25,etiqueta)
        end
    end
end
xlim([-4,makespan(3)+5])
grafica=gca;
set(grafica,'YTick',[])
set(grafica, 'Ydir', 'reverse')
x=[0,0];
y=[0,numeroMaquinas*h];
hold on
plot(x,y,'k','LineWidth',1.25)
hold off
end

%Compare two fuzzy numbers
function [mayor,pos] = mayor_difuso(z1,z2)
mayor=z1;
pos=1;
if sum(z1==z2)==3
    return
end
%Criterio 1, mayor centroide
c1z1=(z1(1)+2*z1(2)+z1(3))/4;
c1z2=(z2(1)+2*z2(2)+z2(3))/4;
if c1z2 > c1z1
    mayor=z2;
    pos=2;
    return
elseif c1z2 == c1z1
    %Criterio 2, mayor valor mas probable
    c2z1=z1(2);
    c2z2=z2(2);
    if c2z2 > c2z1
        mayor=z2;
        pos=2;
        return
    elseif c2z2 == c2z1
        %Criterio 3, mayor diferencia entre cotas
        c3z1=z1(3)-z1(1);
        c3z2=z2(3)-z2(1);
        if c3z2 > c3z1
            mayor=z2;
            pos=2;
            return
        end
    end
end
end

%Sorted indices of fuzzy values
function [indices]=ordenarValoresDifusos(PobValores)
Valores=PobValores;
num_val=size(PobValores,1);
indices=1:num_val;
while(true)
    num_cambios=0;
    for i=1:num_val-1
        z1=Valores(i,:);
        z2=Valores(i+1,:);
        [~,pos] = mayor_difuso(z1,z2);
        if pos==1 && sum(z1-z2==0)<3
            Valores(i,:)=z2;
            Valores(i+1,:)=z1;
            aux=indices(i);
            indices(i)=indices(i+1);
            indices(i+1)=aux;
            num_cambios=num_cambios+1;
        end
    end
    if num_cambios==0
        break
    end
end
end

