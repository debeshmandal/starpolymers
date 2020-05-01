# free energy plotter

def read_data(path):
        data = pd.read_csv(path, skiprows=4, header=None,
                           delim_whitespace=True).rename(
                               columns={
                                   0:'bin',
                                   1:'xi',
                                   2:'count',
                                   3:'p'})

class FreeEnergy():
    
    def __init__(self, ID):
        self.ID = ID
        self.path = 'results/{}/xi.hist'.format(ID)

    def change_path(self, new_path):
        self.path = new_path

    def probability(self):
        data = read_data(self.path)
        return data

    def calculate_energy(self, kBT=1.2):
        data = self.read_data(self)
        p = data['p']
        F = -kBT * np.array([math.log10(num) for num in p])
        data['F']=F
        return data

    def plot_energy(self, wind):
        data = self.calculate_energy(self)
        roller = data.rolling(wind, center=True)
        fig, ax = plt.subplots(1)
        ax.plot(data['xi'], r.mean()['F'], c='red')
        ax.bar(data['xi'], data['F'], alpha=0.4, width=0.01)
        ax.fill_between(data['xi'], r.mean()['F']+r.std()['F'],
                r.mean()['F']-r.std()['F'], alpha=0.4, color='red')
        plt.xlabel(r'$\xi$ [$\sigma$]')
        plt.ylabel(r'F($\xi$) [$k_BT$]')
