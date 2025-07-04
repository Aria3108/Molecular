import numpy as np
from collections import defaultdict
import math

class HMMProfileAnalyzer:
    def __init__(self, sequences):
        """
        Inicializa el analizador HMM Profile con secuencias alineadas
        
        Args:
            sequences: Lista de cadenas alineadas
        """
        self.sequences = sequences
        self.num_sequences = len(sequences)
        self.alignment_length = len(sequences[0]) if sequences else 0
        
        # Aminoácidos estándar
        self.amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                           'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        
        # Probabilidades de background equiprobables
        self.background_prob = {aa: 1/20 for aa in self.amino_acids}
        
        # Variables para almacenar resultados
        self.match_columns = []
        self.num_states = 0
        self.emission_probs = {}
        self.transition_probs = {}
        
    def identify_match_columns(self):
        """
        Identifica columnas que serán estados match (más de la mitad son residuos, no gaps)
        """
        match_columns = []
        
        for col_idx in range(self.alignment_length):
            residue_count = 0
            for seq in self.sequences:
                if col_idx < len(seq) and seq[col_idx] != '-':
                    residue_count += 1
            
            # Si más de la mitad son residuos (no gaps)
            if residue_count > self.num_sequences / 2:
                match_columns.append(col_idx)
        
        self.match_columns = match_columns
        self.num_states = len(match_columns)
        
        print(f"Columnas identificadas como estados match: {[i+1 for i in match_columns]}")
        print(f"Número de estados match: {self.num_states}")
        
        return match_columns
    
    def calculate_emission_probabilities(self):
        """
        Calcula las probabilidades de emisión para cada estado match
        Maneja frecuencia cero usando add-one smoothing
        """
        emission_probs = {}
        
        for state_idx, col_idx in enumerate(self.match_columns):
            state_name = f"M{state_idx + 1}"
            emission_probs[state_name] = {}
            
            # Contar frecuencias de cada aminoácido en la columna
            aa_counts = defaultdict(int)
            total_residues = 0
            
            for seq in self.sequences:
                if col_idx < len(seq) and seq[col_idx] != '-':
                    aa = seq[col_idx]
                    if aa in self.amino_acids:
                        aa_counts[aa] += 1
                        total_residues += 1
            
            # Aplicar add-one smoothing para manejar frecuencia cero
            for aa in self.amino_acids:
                count = aa_counts[aa]
                # Fórmula: (count + 1) / (total + 20)
                prob = (count + 1) / (total_residues + 20)
                emission_probs[state_name][aa] = prob
            
            # Mostrar probabilidades para el estado
            print(f"\nProbabilidades de emisión para {state_name}:")
            for aa in self.amino_acids:
                if emission_probs[state_name][aa] > 0.05:  # Mostrar solo prob significativas
                    print(f"  {aa}: {emission_probs[state_name][aa]:.4f}")
        
        self.emission_probs = emission_probs
        return emission_probs
    
    def calculate_transition_probabilities(self):
        """
        Calcula las probabilidades de transición entre estados
        Considera transiciones Match->Match, Match->Delete, Match->Insert
        """
        transition_probs = {}
        
        for state_idx in range(self.num_states):
            state_name = f"M{state_idx + 1}"
            transition_probs[state_name] = {}
            
            # Contadores para diferentes tipos de transiciones
            match_to_match = 0
            match_to_delete = 0
            match_to_insert = 0
            
            # Analizar transiciones desde este estado
            if state_idx < self.num_states - 1:  # No es el último estado
                current_col = self.match_columns[state_idx]
                next_col = self.match_columns[state_idx + 1]
                
                for seq in self.sequences:
                    # Si la secuencia tiene residuo en la columna actual
                    if (current_col < len(seq) and seq[current_col] != '-'):
                        
                        # Verificar qué pasa en la siguiente columna match
                        if next_col < len(seq) and seq[next_col] != '-':
                            # Hay residuo en la siguiente columna match
                            # Verificar si hay inserciones entre medio
                            insertion_found = False
                            for intermediate_col in range(current_col + 1, next_col):
                                if (intermediate_col < len(seq) and 
                                    seq[intermediate_col] != '-' and 
                                    intermediate_col not in self.match_columns):
                                    insertion_found = True
                                    break
                            
                            if insertion_found:
                                match_to_insert += 1
                            else:
                                match_to_match += 1
                        else:
                            # No hay residuo en la siguiente columna match
                            match_to_delete += 1
            
            # Calcular probabilidades con add-one smoothing
            total_transitions = match_to_match + match_to_delete + match_to_insert
            
            if total_transitions > 0:
                transition_probs[state_name][f"M{state_idx + 2}"] = (match_to_match + 1) / (total_transitions + 3)
                transition_probs[state_name][f"D{state_idx + 2}"] = (match_to_delete + 1) / (total_transitions + 3)
                transition_probs[state_name][f"I{state_idx + 1}"] = (match_to_insert + 1) / (total_transitions + 3)
            else:
                # Si no hay transiciones, usar probabilidades uniformes
                transition_probs[state_name][f"M{state_idx + 2}"] = 1/3
                transition_probs[state_name][f"D{state_idx + 2}"] = 1/3
                transition_probs[state_name][f"I{state_idx + 1}"] = 1/3
        
        self.transition_probs = transition_probs
        
        # Mostrar probabilidades de transición
        print("\nProbabilidades de transición:")
        for state, transitions in transition_probs.items():
            print(f"{state}:")
            for target, prob in transitions.items():
                print(f"  -> {target}: {prob:.4f}")
        
        return transition_probs
    
    def build_hmm_profile(self):
        """
        Construye el profile HMM completo
        """
        print("=== CONSTRUCCIÓN DEL PROFILE HMM ===")
        print(f"Secuencias de entrada:")
        for i, seq in enumerate(self.sequences):
            print(f"  {i+1}: {seq}")
        
        print(f"\nLongitud del alineamiento: {self.alignment_length}")
        print(f"Número de secuencias: {self.num_sequences}")
        
        # Paso 1: Identificar columnas match
        self.identify_match_columns()
        
        # Paso 2: Calcular probabilidades de emisión
        self.calculate_emission_probabilities()
        
        # Paso 3: Calcular probabilidades de transición
        self.calculate_transition_probabilities()
        
        return {
            'num_states': self.num_states,
            'match_columns': self.match_columns,
            'emission_probs': self.emission_probs,
            'transition_probs': self.transition_probs
        }
    
    def analyze_sequence(self, test_sequence):
        """
        Analiza si una secuencia de prueba pertenece a la familia
        usando el profile HMM construido
        """
        if not self.emission_probs:
            print("Error: Debe construir el profile HMM primero")
            return None
        
        # Calcular score usando el método de ventana deslizante
        max_score = float('-inf')
        best_position = -1
        
        window_size = self.num_states
        
        for pos in range(len(test_sequence) - window_size + 1):
            window = test_sequence[pos:pos + window_size]
            score = 0
            
            for i, aa in enumerate(window):
                state_name = f"M{i + 1}"
                if state_name in self.emission_probs and aa in self.emission_probs[state_name]:
                    # Score = log(P_ij / b_i) donde P_ij es prob de emisión y b_i es background
                    emission_prob = self.emission_probs[state_name][aa]
                    background_prob = self.background_prob.get(aa, 1/20)
                    score += math.log(emission_prob / background_prob)
            
            if score > max_score:
                max_score = score
                best_position = pos
        
        print(f"\nAnálisis de secuencia: {test_sequence}")
        print(f"Mejor posición: {best_position}")
        print(f"Score máximo: {max_score:.4f}")
        
        return {
            'sequence': test_sequence,
            'best_position': best_position,
            'max_score': max_score
        }

# Función principal para ejecutar el análisis
def main():
    # Secuencias del ejemplo de clase
    sequences = [
        "VGA--HAGEY",
        "V----NVDEV", 
        "VEA--DVAGH",
        "VKG------D",
        "VYS--TYETS",
        "FNA--NIPKH",
        "IAGADNGAGY"
    ]
    
    # Crear analizador
    analyzer = HMMProfileAnalyzer(sequences)
    
    # Construir el profile HMM
    hmm_profile = analyzer.build_hmm_profile()
    
    # Mostrar resumen
    print("\n=== RESUMEN DEL PROFILE HMM ===")
    print(f"Número de estados match: {hmm_profile['num_states']}")
    print(f"Columnas consideradas como estados match: {[i+1 for i in hmm_profile['match_columns']]}")
    
    # Ejemplo de análisis de una secuencia de prueba
    test_sequence = "VGAHAGEY"
    analyzer.analyze_sequence(test_sequence)
    
    return analyzer

if __name__ == "__main__":
    analyzer = main()