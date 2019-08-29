from astropy.io import fits #pacote para trabalhar com imagens em formato fits
import numpy as np #pacote para trabalhar com matrizes
import matplotlib.pyplot as plt
#
def mag_galaxia():
    '''Função a ser execturada'''
#
    centro = [1492, 1312]
#
    img = processar_imagens('img.fits', 'mbias.fits', 'flat.fits')
    tamanho = [img[0].shape[0], img[0].shape[1]]
    raio_gal = 20
    raio_ceu = 45
#
    #vai criar uma matriz em que todos os valores fora da região ao redor da galáxia vai ser nulos
    mascara_gal = create_circular_mask(tamanho[0], tamanho[1], centro, raio_gal)
    galaxia = img[0]*mascara_gal
#
    #agora, a intenção é criar uma região ao redor da galáixa em forma de anel que exclua a vizinhança próxima da galáxia
    mascara_externo = create_circular_mask(tamanho[0], tamanho[1], centro, raio_ceu)
    interno = mascara_externo*img[0]
    mascara_interna = create_circular_mask_interno(tamanho[0], tamanho[1], centro, raio_gal+15)
    ceu = mascara_interna*interno
#
    mascara_bias = create_circular_mask(tamanho[0], tamanho[1], centro, raio_gal)
    bias = img[1][1]*mascara_bias
#
#
    mediana_ceu = np.median(ceu[ceu != 0]) #calcula a media dentro do anel estipulado para o céu
#
    img_final = galaxia - mediana_ceu; img_final[img_final < 0] = 0
#
    contagem_gal = np.sum(img_final) #calcula a contagem dentro do raio estipulado para a galáxia
    contagem_bias = np.sum(bias[bias != 0])
    contagem_ceu = np.sum(ceu)
    contagem_estrela_padrao = mag_padrao(img) #calcula a contagem da estrela padrão
#
    magnitude = 18 + 2.5*np.log(contagem_estrela_padrao/contagem_gal) #calcula a magnitude da galáxia tendo a estrela padrão como mag de referência
#
    sinal_ruido = contagem_gal/np.sqrt((np.sum(galaxia[galaxia != 0]))+contagem_bias+contagem_ceu) #calcula o S/N da galáxia
#
    erro_mag = np.abs(2.5*np.log10(1 + 1/sinal_ruido))
#
    print("A magnitude da galáxia é: ", magnitude)
    print("O erro da magnitude é: ", erro_mag)
    print("O S/N da galáxia é: ", sinal_ruido)

##================================================================##===========================================================##
##================================================================##===========================================================##
##================================================================##===========================================================##


def processar_imagens(caminho_imagem, caminho_bias, caminho_flat):
    IMG = []
    HDR = []
    caminhos = [caminho_imagem, caminho_bias, caminho_flat]
    for i in range(len(caminhos)):
        img, hdr = fits.getdata(caminhos[i], header=True) #salva em img apenas a matriz de dados
        img = np.array(img, dtype='Float64')   # converte os dados para float64
        IMG.append(img) #adiciona os dados atuais na lista geral de dados
        HDR.append(hdr)
    img_final = (IMG[0] - IMG[1])/IMG[2]
    return img_final, IMG, HDR

##================================================================##===========================================================##
##================================================================##===========================================================##
##================================================================##===========================================================##


def mag_padrao(img):
    tamanho = [img[0].shape[0], img[0].shape[1]]
#
    raio_ceu = 45
    raio_estrela_padrao = 13.5
    raio_ceu_padrao = 30
    centro = [1492, 1312]
#
    #vai criar uma matriz em que todos os valores fora da região ao redor da estrela vai ser nulos
    mascara_estrela_padrao = create_circular_mask(tamanho[0], tamanho[1], centro, raio_estrela_padrao)
    estrela_padrao = img[0]*mascara_estrela_padrao
#
    #agora, a intenção é criar uma região ao redor da estrela em forma de anel que exclua a vizinhança próxima da estrela
    mascara_externo = create_circular_mask(tamanho[0], tamanho[1], centro, raio_ceu)
    interno = mascara_externo*img[0]
#
    mascara_interna = create_circular_mask_interno(tamanho[0], tamanho[1], centro, raio_estrela_padrao+7.5)
    ceu_padrao = mascara_interna*interno
#
    mediana_ceu_padrao = np.median(ceu_padrao[ceu_padrao != 0]) #calcula a media dentro do anel estipulado para o céu
    img_final_padrao = estrela_padrao - mediana_ceu_padrao
    contagem_estrela_padrao = np.sum(estrela_padrao) #calcula a contagem dentro do raio estipulado para a estrela
    return contagem_estrela_padrao

##================================================================##===========================================================##
##================================================================##===========================================================##
##================================================================##===========================================================##


def create_circular_mask(h, w, center, radius):
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center <= radius
    return mask

def create_circular_mask_interno(h, w, center, radius):
    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)
    mask = dist_from_center >= radius
    return mask
##================================================================##===========================================================##
