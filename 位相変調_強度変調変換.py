# 200THzを30GHzで変調する

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rayleigh

class Modulador:

    def __init__(self, Fc, Fs, Tb):
        self.Fc = Fc
        self.Fs = Fs
        self.Ts = 1/Fs
        self.Tb = Tb  # Largura de cada símbolo.
        self.L = int(Tb * Fs)  # Número de amostras para o tempo de cada símbolo.

    def processar(self, simbolos):
        # Gera a onda quadrada (cada Tb dura L amostras).
        ondaq = np.repeat(simbolos, self.L)

        # Gera a portadora.
        t = np.linspace(0, ondaq.size, ondaq.size)
        portadora = np.cos(2 * np.pi * self.Fc * t)
        sinalm = np.multiply(ondaq, portadora)

        return (ondaq, sinalm)


class Demodulador:

    def __init__(self, modulador):
        self.modulador = modulador

    def processar(self, sinalm):

        # Estágio 1

        # Gera a portadora.
        t = np.linspace(0, sinalm.size, sinalm.size)
        portadora = np.cos(2 * np.pi * self.modulador.Fc * t)
        sinald = np.multiply(sinalm, portadora)

        # Integra para melhorar o formato da onda (opcional).
        sinali = np.convolve(sinald, np.ones(self.modulador.L))

        # Remove o atraso de self.modulador.L - 1 amostras.
        sinali = sinali[int(self.modulador.L) - 1::]

        # Estágio 2 (decisor)

        # Decide se é 1 ou -1 baseado no limiar 0.
        positivos = (sinali > 0)
        negativos = np.logical_not(positivos)

        ondaq = np.empty(sinali.size)
        ondaq[positivos] = 1
        ondaq[negativos] = -1

        # Faz uma subamostragem para obter os símbolos.
        simbolos = ondaq[::self.modulador.L]

        return (sinald, sinali, ondaq, simbolos)


class Canal:

    def __init__(self, SNR, taps, Fd):
        self.SNR = SNR
        self.taps = taps
        self.Fd = Fd

    def processar(self, sinal):


        # Aplicando ruído gaussiano branco.
        potencia_sinal = np.sum(np.square(sinal))/sinal.size

        # Gerando ruído gaussiano branco (média = 0, variancia = potencia do awgn).
        potencia_ruido = potencia_sinal / self.SNR
        desvio_padrao = np.sqrt(potencia_ruido)
        ruido_gaussiano = np.random.normal(0, desvio_padrao, sinal.size)

        # Aplica o ruído ao sinal.
        sinal_ruidoso = sinal + ruido_gaussiano

        # Atenua o sinal.
        return sinal_ruidoso

# Configurações.

N = 100  # Número de símbolos a ser enviado.
Fc =30 * 10 ** 9  # Frequência da portadora.
Fs = 4 * Fc  # Frequência de amostragem.
Tb = 0.1  # Largura de cada símbolo (em seg).

SNRdB = 100  # Potência do sinal é duas vezes o dobro da potência do ruído.
SNR = 10.0 ** (SNRdB/10.0)

Fd = 5  # Frequência de Doppler para o canal de múltiplos caminhos.

TAPS = 5  # Número de elementos (caminhos) do canal.

# Cria o modulador e o demodulador.
modulador = Modulador(Fc, Fs, Tb)
demodulador = Demodulador(modulador)
canal = Canal(SNR, TAPS, Fd)

# Dados a serem enviados.
dados = np.random.choice(np.array([0, 1]), size=(N))

# Criando símbols para o BPSK (-1 e 1).
simbolos_enviados = 2*dados - 1

# Modula os símbolos.
(ondaq_enviada, sinalm) = modulador.processar(simbolos_enviados)

# Processa pelo canal.
sinalc = canal.processar(sinalm)

# Demodula o sinal recebido.
(sinald,  sinali, ondaq_recebida, simbolos_recebidos) = demodulador.processar(sinalc)

dados_recebidos = ((simbolos_recebidos + 1)/2).astype(int)

# Calculando os erros de decisão.
num_erros = np.sum(simbolos_enviados != simbolos_recebidos)
BER = num_erros/simbolos_enviados.size

print('Do total de {} bits, {} foram decodificados de formada errada.'.format(
    simbolos_enviados.size, num_erros
))
print('BER: {}'.format(BER))

# Exibindo gráficos do transmissor.
f1, (f1_ax1, f1_ax2, f1_ax3) = plt.subplots(3)
f1.suptitle('Sinal enviado a partir do transmissor', fontsize=14)
f1_ax1.stem(dados)
f1_ax1.set_title('Bits enviados')
f1_ax2.plot(ondaq_enviada)
f1_ax2.set_title('Onda quadrada gerada a partir dos símbolos')
f1_ax3.plot(sinalm)
f1_ax3.set_title('Sinal modulado')
f1.subplots_adjust(hspace=1)

# Exibindo gráficos do receptor.
f2, (f2_ax1, f2_ax2, f2_ax3, f2_ax4, f2_ax5) = plt.subplots(5)
f2.suptitle('Sinal recebido no receptor.', fontsize=14)
f2_ax1.plot(sinalc)
f2_ax1.set_title('Sinal recebido do canal')
f2_ax2.plot(sinald)
f2_ax2.set_title('Sinal demodulado')
f2_ax3.plot(sinali)
f2_ax3.set_title('Sinal após integração')
f2_ax4.plot(ondaq_recebida)
f2_ax4.set_title('Onda quadrada recebida')
f2_ax5.stem(dados_recebidos)
f2_ax5.set_title('Dados recebidos')
f2.subplots_adjust(hspace=1)
plt.show()