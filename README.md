# An Introduction to Economic Dynamics - Implementation

This repository contains Julia code that implements the concepts discussed in the book **"An Introduction to Economic Dynamics"** by **Srinivas Raghavendra and Petri T. Piiroinen**.

## ⚠️ Disclaimer

This repository is an independent project created for educational purposes. It is not associated with the authors of the book.

## 📁 Project Structure

```
repo/
├── env/          # Code for initializing the local environment
├── figures/      # Visualizations generated by the code
├── src/          # Source code implementing economic dynamics
├── .gitignore
├── Manifest.toml
├── Project.toml
└── README.md     # Project documentation
```

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/ManullangJihan/Economic-Dynamics-Modeling.git
```

### 2. Initialize Local Environment

The local environment is set up using the `activate_env.jl` file inside the `env` folder. Each script in the `src` directory begins with the following line:

```julia
# Activate environment and install dependencies defined in Project.toml
include("../env/activate_env.jl")
```

### 3. Running Code

* The main source code is in the `src` folder.
* The visualizations will be saved in the `figures` folder.

### 📚 Reference

**"An Introduction to Economic Dynamics"** by **Srinivas Raghavendra and Petri T. Piiroinen**.

