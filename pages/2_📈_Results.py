import streamlit as st
from data_manager import load_data
from sidebar import sidebar

st.title("📈 Results")
sidebar("results_page")

# Читаем финальное состояние файла
final_data = load_data()

st.write("### Data from your current session file:")
st.json(final_data)

# Тут вызываешь свой алгоритм, передавая ему final_data