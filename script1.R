#setwd("C:\\Users\\sonia\\source\\repos\\RProject1\\RProject1")
library(signal, warn.conflicts = FALSE)
library(seewave)

#Leemos el fichero
x <- read.csv("samples.csv", header = TRUE);

#Array con las amplitudes (voltios)
x1 = x$volt;#[1:500];
#N?mero de muestras
N = length(x1);
#Frecuencia de muestreo
fs = N / 10;
#Array con los tiempos
t = x$time;#[1:500];


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
x6 = x6[14 + (1:N)];
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
poss_reg = (x6 > ((thresh * max_h)));

#Array con la lista de unos (las esquinas por la izquierda)
left = which(diff(c(0, poss_reg)) %in% 1);
#Array con la lista de -1 (las esquinas por la derecha)
right = which(diff(c(poss_reg, 0)) %in% -1);

left = left - (delay);
right = right - (delay);




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
    Q_value = min(x1[left[i]:R_loc]);
    indexOfQ = which(x1[left[i]:right[i]] == Q_value);
    Q_time = t[left[i] + indexOfQ[1]];
    Q_loc = left[i] + (Q_value - 1);

    #PuntoS
    S_value = min(x1[left[i]:right[i]]);
    indexOfS = which(x1[left[i]:right[i]] == S_value);
    S_time = t[left[i] + indexOfS[1]];
    S_loc = left[i] + (S_value - 1);

    #Almacenamos el momento del tiempo en el que ocurre cada evento
    timeValues[i,] <- c(R_time, Q_time, S_time);

    #Almacenamos los valores QRS
    qrsValues[i,] <- c(R_value, Q_value, S_value);


}


#Array que contiene los colores de los puntos
colors <- c("green", "blue", "red")


plot(t, x1, type = "l", xlab = "time", ylab = "volts")
legend(0, 0.8, c("Q", "R", "S"), pch = c(15, 15, 15), lwd = c(2.5, 2.5, 2.5), col = c("blue", "green", "red"))
for (i in 1:length(left)) {
    for (x in 1:3) {
        points(timeValues[i, x], qrsValues[i, x], pch = 15, col = colors[x]);
    }
}

