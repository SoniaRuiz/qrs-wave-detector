setwd("C:\\Users\\sonia\\source\\repos\\RProject1\\RProject1")
library(signal, warn.conflicts = FALSE)
library(seewave)

#Leemos el fichero
x <- read.csv("samples7.csv", header = TRUE);

#Array con las amplitudes (voltios)
x1 = x$volt;
#Número de muestras
N = length(x1);
#Frecuencia de muestreo
fs = N / 10;
#Array con los tiempos
t = x$time;


x1 = x1 - mean(x1);
#x1 = x1 / max(abs(x1));

#Filtro de paso bajo (Low Pass Filtering)
#LPF(1 - z ^ -6) ^ 2 / (1 - z ^ -1) ^ 2
b = c(1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1);
a = c(1, -2, 1);
x2 = filter(b, a, x1);
x2 = x2 / max(abs(x2));
delay = 6;



# High Pass filter H(z) = (-1 + 32 z ^ (-16) + z ^ (-32)) / (1 + z ^ (-1))
b = c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1);
a = c(1, 1);
x3 = filter(b, a, x2);
x3 = x3 / max(abs(x3));
delay = delay + 16;


#derivative filter H(z) = (1 / 8 T)(-z ^ (-2) - 2 z ^ (-1) + 2 z + z ^ (2))
b = c(1, 2, 0, -2, -1) * (1 / 8) * fs;
x4 = conv(b, x3);
x4 = x4 / max(x4);


# Squaring nonlinearly enhance the dominant peaks
x5 = x4 ^ 2;


#Movemos la ventana de integracion (Moving Window Integration)
h = matrix(1, 1, 31) / 31;
x6 = conv(x5, h);
x6 = x6[15 + (1:N)];
x6 = x6 / max(abs(x6));

#x6 = conv(x5, array(1, dim = c(1, round(0.150 * fs))) / round(0.150 * fs));
#delay = delay + round(0.150 * fs) / 2;
#plot((1:length(x6) - 1) / fs, x6)


#Busqueda del complejo QRS (Find QRS Points Which it is different than Pan - Tompkins algorithm)
#Maximo
max_h = max(x6);
#Media
thresh = mean(x6);
#TRUE cuando detecta valores que contienen el complejo QRS
poss_reg = (x6 > ((thresh * max_h)*2));

#Array con la lista de unos (las esquinas por la izquierda)
left = which(diff(c(0, poss_reg)) %in% 1);
#Array con la lista de -1 (las esquinas por la derecha)
right = which(diff(c(poss_reg, 0)) %in% -1);

left = left - (delay);
right = right - (delay);

for (i in 1:length(left)) {
    if (left[i] < 0)
        left[i] = 0;
    }


#Matriz con los puntos QRS
qrsValues <- matrix(0, nrow = length(left), ncol = 3)
#Matriz con los tiempos
timeValues <- matrix(0, nrow = length(left), ncol = 3)



for (i in 1:length(left)) {

    #PuntoR
    R_value = max(x1[left[i]:right[i]]);
    indexOfR = which(x1[left[i]:right[i]] == R_value);
    R_time = t[left[i] + indexOfR[1]];
    R_loc = left[i] + (indexOfR[1] - 1);

    #PuntoQ
    Q_value = min(x1[left[i]:R_loc])[1];
    indexOfQ = which(x1[left[i]:R_loc] == Q_value);
    Q_time = t[left[i] + indexOfQ[1]];
    Q_loc = left[i] + (Q_value - 1);

    #PuntoS
    S_value = min(x1[R_loc:right[i]]);
    indexOfS = which(x1[left[i]:right[i]] == S_value);
    if (length(indexOfS) > 1) {
        S_time = t[left[i] + indexOfS[length(indexOfS)]];
    }
    else {
        S_time = t[left[i] + indexOfS];
    }
    S_loc = left[i] + (S_value - 1);

    #Almacenamos el momento del tiempo en el que ocurre cada evento
    timeValues[i,] <- c(R_time, Q_time, S_time);

    #Almacenamos los valores QRS
    qrsValues[i,] <- c(R_value, Q_value, S_value);

}
#Array que contiene los colores de los puntos
colors <- c("green", "blue", "red")

plot(t,x1, type = "l", xlab = "time", ylab = "volts", ylim = c(-0.3, 1.2))

for (i in 1:length(left)) {
    for (x in 1:3) {
        #points(timeValues[i, x], qrsValues[i, x], pch = 15, col = colors[x]);
    }
}


##################################################
##################################################
if (N >= 3000) {
    totalBeatsPerMinute = 0;
    for (y in 1:length(left)) {
        if (y != length(left)) {
            interval = timeValues[y + 1, 1] - timeValues[y, 1];
            beatsPerMinute = (1 / interval) * 60;
            totalBeatsPerMinute = totalBeatsPerMinute + beatsPerMinute;
        }
    }
    totalBeatsPerMinute = mean(totalBeatsPerMinute / 10);
    text(x = 0.5, y = 1.2, labels = paste("Ppm:", round(totalBeatsPerMinute, digits = 2), sep = " "));
}
##################################################
##################################################

#Matriz con los puntos P
pValues <- matrix(0, nrow = length(left) + 1, ncol = 3)
#Matriz con los tiempos de P
pTimeValues <- matrix(0, nrow = length(left) + 1, ncol = 3)

#Matriz con los puntos T
tValues <- matrix(0, nrow = length(left) + 1, ncol = 3)
#Matriz con los tiempos de T
tTimeValues <- matrix(0, nrow = length(left) + 1, ncol = 3)

offset = 0;
u = 1;
for (i in 1:(length(left) + 1)) {
    #Si estoy en la primera posicion
    if (i == 1) {
        g1 <- x1[1:left[i]];
        offset = 0;
    }
    else {
        offset = right[i - 1];
        #Si estoy en la última posición
        if (i == length(left) + 1) {
            g1 <- x1[right[i - 1]:length(x1)];
        }
        else
            g1 <- x1[right[i - 1]:left[i]];
        }

    g1 = g1 / max(abs(g1));
    g1 = g1 - mean(g1);

    #Filtro pasa-banda
    g2 <- seewave::bwfilter(g1, fs, from = 0.5, to = 10, bandpass = TRUE, output = "matrix");
    g2 <- rev(seewave::bwfilter(rev(g2), fs, from = 0.5, to = 10, bandpass = TRUE, output = "matrix"));

    ##########################
    ## CÁLCULO DE LÍMITES SEGÚN PPM
    #########################
    lowLimit = (totalBeatsPerMinute * 55) / 60;
    highLimit = (totalBeatsPerMinute * 110) / 60;

    #Paso 3: Media móvil
    if (length(g2) < lowLimit)
        MAPeak <- stats::filter(g2, rep(1 / lowLimit, length(g2)), method = c("convolution"), sides = 2, circular = TRUE)
    else
        MAPeak <- stats::filter(g2, rep(1 / lowLimit, lowLimit), method = c("convolution"), sides = 2, circular = TRUE)

    if (length(g2) < highLimit)
        MAPwave <- stats::filter(g2, rep(1 / highLimit, length(g2)), method = c("convolution"), sides = 2, circular = TRUE)
    else
        MAPwave <- stats::filter(g2, rep(1 / highLimit, highLimit), method = c("convolution"), sides = 2, circular = TRUE)


    #Paso 4: Seleccionamos los bloques de interés
    #Cuando la amplitud de la primera media móvil (MAPeak) es mayor que la amplitud de la segunda media móvil (MAPwave),
    #se selecciona como bloque de interés
    blocks <- array(-1, dim = c(1, length(MAPeak)))
    for (y in 1:length(MAPeak)) {
        if ((MAPeak[y] > MAPwave[y])) {
            blocks[y] <- 0.25;
        }
        else
            blocks[y] <- 0;
        }

    #Paso 5: Eliminamos los bloques de ruido
    numBlocksGreaterThanZero = blocks > 0;
    leftBlocks = which(diff(c(0, numBlocksGreaterThanZero)) %in% 1);
    rightBlocks = which(diff(c(numBlocksGreaterThanZero, 0)) %in% -1);
    if (length(leftBlocks) == 0 || length(rightBlocks) == 0)
        break;
    for (y in 1:length(leftBlocks)) {
        if (length(blocks[leftBlocks[y]:rightBlocks[y]]) < (lowLimit * 0.5)) {
            blocks[leftBlocks[y]:rightBlocks[y]] <- 0;
        }
    }



    #Paso 6. Para determinar si cada bloque detectado en la sección Ri-Ri+1 contiene una onda P o T, se cuenta el número de bloques de ese intervalo.
    numBlocksGreaterThanZero = blocks > 0;
    leftBlocks = which(diff(c(0, numBlocksGreaterThanZero)) %in% 1);
    rightBlocks = which(diff(c(numBlocksGreaterThanZero, 0)) %in% -1);
    #A continuación se extrae el máximo para distinguir la onda P de la T y del ruido
    #Si se ha detectado 0 bloques: error (salimos del bucle
    if (length(leftBlocks) == 0 || length(rightBlocks) == 0)
        break;

    #Si se ha detectado 1 o más bloques:
    ### POSICIONES DE LAS R
    #Si estoy en el primer bloque
    if (i == 1) {
        RTimeValue1 = 0.00;
        RTimeValue2 = timeValues[i, 1];
    }
    #Si estoy en el último bloque (ya no hay más complejos QRS)
    else if (i == length(left) + 1) {
        RTimeValue1 = timeValues[i - 1, 1];
        RTimeValue2 = t[length(t)];
    }
    #Para el resto de bloques
    else {
        RTimeValue1 = timeValues[i - 1, 1];
        RTimeValue2 = timeValues[i, 1];
    }


    max = 0;
    for (y in 1:length(leftBlocks)) {

        ## PUNTO POR LA IZQUIERDA Y POR LA DERECHA DE CADA BLOQUE DE INTERÉS
        leftBlock = offset + leftBlocks[y];
        rightBlock = (offset + rightBlocks[y]) - 1;


        ## MÁXIMO VALOR DE ESE BLOQUE DE INTERÉS
        maxBlockInLoop = max(x1[leftBlock:rightBlock]);



        ## POSICIÓN DEL VALOR MÁXIMO DENTRO DE LA MATRIZ ORIGINAL
        #Si estoy en la primera región Ri-Ri+1
        if (i == 1) {
            indexOfMax = (which(x1[1:rightBlock] == maxBlockInLoop));
        }
        #Si estoy en la última región Ri-Ri+1
        else if (i == length(left) + 1) {

            indexOfMax = offset + (which(x1[right[i - 1]:rightBlock] == maxBlockInLoop)) - 1;

        }
        #Si estoy en cualquier otra región
        else {
            indexOfMax = offset + (which(x1[right[i - 1]:left[i]] == maxBlockInLoop));
        }



        if (maxBlockInLoop > max) {
            max = maxBlockInLoop;
            timeValueMax = t[indexOfMax[1]];
            distanceToR = timeValueMax - RTimeValue1;

            ## VALORES MÁXIMOS DE LAS DISTANCIAS A LAS R
            RTmin = (RTimeValue2 - RTimeValue1) * 0.111;
            RTmax = (RTimeValue2 - RTimeValue1) * 0.583;

            if (distanceToR > RTmin && distanceToR < RTmax) {
                #Almacenamos los valores T
                tValues[i,] <- c(maxBlockInLoop, x1[leftBlock], x1[rightBlock]);
                #Almacenamos los tiempos de P
                tTimeValues[i,] <- c(timeValueMax, t[leftBlock], t[rightBlock]);
            }
        }
        else {

            #Solo analizo si no estoy en la último posición. Si estoy en la última posición es que ya no hay más complejos QRS, 
            #por lo que no puedo detectar la onda P con seguridad
            if (i < length(left) + 1) {
                timeValueMax = t[indexOfMax[length(indexOfMax)]];
                distanceToR = RTimeValue2 - timeValueMax;

                if ((RTimeValue2 - t[rightBlock]) > (0.055 * 0.5) && (RTimeValue2 - t[rightBlock]) < (0.47 * 0.5)) {
                    #Almacenamos los valores P
                    pValues[i,] <- c(maxBlockInLoop, x1[leftBlock], x1[rightBlock]);
                    #Almacenamos los tiempos de P
                    pTimeValues[i,] <- c(timeValueMax, t[leftBlock], t[rightBlock]);
                }
            }
        }
    }
}


for (i in 1:(length(left) + 1)) {
    #print(i)
    for (x in 1:3) {
        if (tTimeValues[i, x] != 0 && tValues[i, x] != 0)
            points(tTimeValues[i, x], tValues[i, x], pch = 15, col = "purple");
        if (pTimeValues[i, x] != 0 && pValues[i, x] != 0)
            points(pTimeValues[i, x], pValues[i, x], pch = 15, col = "orange");
        }
}
#legend(0, 1.2, c("Q", "R", "S"), lty = c(1, 1, 1), lwd = c(2.5, 2.5, 2.5), col = c("black", "orange", "purple"))
legend(0, 1, c("QRS Removal"), lty = c(1), lwd = c(2.5), col = c("red"))
#lines((((1:length(blocks) - 1) + offset) / fs), blocks, type = "l", col = "red");
