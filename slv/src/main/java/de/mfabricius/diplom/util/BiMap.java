package de.mfabricius.diplom.util;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;


/**
 * This class provides a simple implementation of a bidirectional {@link Map}. It can only store unique mappings. If a
 * new key would point to an existing value, the Map will reject the new key and throw an {@link IllegalStateException}.
 * 
 * @author Martin Fabricius
 *
 * @param <K>
 *            type of the keys
 * @param <V>
 *            type of the values
 */
public class BiMap<K, V> implements Map<K, V> {

    private final HashMap<K, V> map = new HashMap<>();
    private final HashMap<V, K> inverse = new HashMap<>();

    @Override
    public V get(Object key) {
        return map.get(key);
    }

    /**
     * Get the associated key for the given value.
     * 
     * @param value
     *            the value to look up
     * @return the associated key
     */
    public K getInverse(V value) {
        return inverse.get(value);
    }

    @Override
    public V put(K key, V value) {
        final V put = map.put(key, value);
        final K put2 = inverse.put(value, key);
        if (put2 != null) {
            map.put(key, put);
            inverse.put(value, put2);
            throw new IllegalStateException("No unique mapping!");
        }
        if (put != null) inverse.remove(put);
        return put;
    }

    /** @return the unmodifiable map of the key-value pairs. */
    public Map<K, V> getMap() {
        return Collections.unmodifiableMap(map);
    }

    /** @return the unmodifiable reverse map of the value-key pairs. */
    public Map<V, K> getInverseMap() {
        return Collections.unmodifiableMap(inverse);
    }

    @Override
    public void clear() {
        map.clear();
        inverse.clear();
    }

    @Override
    public int size() {
        return map.size();
    }

    @Override
    public boolean isEmpty() {
        return map.isEmpty();
    }

    @Override
    public boolean containsKey(Object key) {
        return map.containsKey(key);
    }

    @Override
    public boolean containsValue(Object value) {
        return inverse.containsKey(value);
    }

    @Override
    public V remove(Object key) {
        final V remove = map.remove(key);
        if (remove != null) inverse.remove(remove);
        return remove;
    }

    @Override
    public void putAll(Map<? extends K, ? extends V> m) {
        m.entrySet().forEach(e -> put(e.getKey(), e.getValue()));
    }

    @Override
    public Set<K> keySet() {
        return map.keySet();
    }

    @Override
    public Collection<V> values() {
        return map.values();
    }

    @Override
    public Set<java.util.Map.Entry<K, V>> entrySet() {
        return map.entrySet();
    }
}